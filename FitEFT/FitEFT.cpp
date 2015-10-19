//
//  FitEFT.cpp
//  FitEFT
//
//  Created by Christopher Jacobsen on 23/09/15.
//  Copyright (c) 2015 Christopher Jacobsen. All rights reserved.
//

#include "FitEFT.h"
#include "RootUtil.h"

// Root includes
#include <TMinuitMinimizer.h>

#include <TFile.h>
#include <TH1.h>
#include <TProfile.h>
#include <TNtupleD.h>

#include <TF1.h>
#include <TGraph.h>
#include <TGraphErrors.h>

#include <TFitResult.h>
#include <TVirtualFitter.h>
#include <TBackCompFitter.h>
#include <Math/Minimizer.h>
#include <Math/MinimizerOptions.h>

#include <HFitInterface.h>
#include <Fit/UnBinData.h>
#include <Foption.h>

#include <Fit/FitUtil.h>
#include <Math/WrappedMultiTF1.h>

using namespace RootUtil;

////////////////////////////////////////////////////////////////////////////////
namespace FitEFT
{

////////////////////////////////////////////////////////////////////////////////

static const ROOT::Math::MinimizerOptions StartupMinimizerOptions;

void RestoreDefaultMinimizerOptions()
{
    using namespace ROOT::Math;

    const MinimizerOptions & opt = StartupMinimizerOptions;

    MinimizerOptions::SetDefaultMinimizer(        opt.MinimizerType()     .c_str(),
                                                  opt.MinimizerAlgorithm().c_str() );
    MinimizerOptions::SetDefaultErrorDef(         opt.ErrorDef()         );
    MinimizerOptions::SetDefaultTolerance(        opt.Tolerance()        );
    MinimizerOptions::SetDefaultPrecision(        opt.Precision()        );
    MinimizerOptions::SetDefaultMaxFunctionCalls( opt.MaxFunctionCalls() );
    MinimizerOptions::SetDefaultMaxIterations(    opt.MaxIterations()    );
    MinimizerOptions::SetDefaultStrategy(         opt.Strategy()         );
    MinimizerOptions::SetDefaultPrintLevel(       opt.PrintLevel()       );
    MinimizerOptions::SetDefaultExtraOptions(     opt.ExtraOptions()     );

    //MinimizerOptions test;
    //test.Print();
}

////////////////////////////////////////////////////////////////////////////////
void ClearLastFitter()
{
    TVirtualFitter * pLastFitter = TVirtualFitter::GetFitter();  // get the global last fitter
    if (pLastFitter)
    {
        // the destructor of TVirtualFitter clears the global last fitter
        // if the object deleted is also the global fitter

        delete pLastFitter;

        pLastFitter = TVirtualFitter::GetFitter();
        if (pLastFitter)
            ThrowError( "Failed to clear last fitter." );
    }
}

////////////////////////////////////////////////////////////////////////////////
TGraph * GraphFromProfile( const TProfile & profile, bool bWithErrors = true )
{
    const Int_t nBins = profile.GetNbinsX();

    TGraphErrors * pGraphErrors =  bWithErrors ? new TGraphErrors( nBins ) : nullptr;
    TGraph *       pGraph       = !bWithErrors ? new TGraph(       nBins ) : pGraphErrors;

    // copy name and titles
    pGraph->SetName(  profile.GetName()  );
    pGraph->SetTitle( profile.GetTitle() );
    pGraph->GetXaxis()->SetTitle( profile.GetXaxis()->GetTitle() );
    pGraph->GetYaxis()->SetTitle( profile.GetYaxis()->GetTitle() );

    Int_t index = 0;
    for (Int_t bin = 1; bin <= nBins; ++bin)
    {
        Double_t nEff = profile.GetBinEffectiveEntries(bin);
        if (nEff <= 0.0)
            continue;

        Double_t xValue = profile.GetBinCenter(bin);
        Double_t yValue = profile.GetBinContent(bin);

        pGraph->SetPoint( index, xValue, yValue );

        if (pGraphErrors)
        {
            Double_t yError = profile.GetBinError(bin);

            pGraphErrors->SetPointError( index, 0.0, yError );    // xError not defined
        }

        ++index;
    }

    pGraph->Set(index);

    return pGraph;
}

////////////////////////////////////////////////////////////////////////////////

struct FitResult
{
    struct Param
    {
        const char *    name            = nullptr;
        double          value           = 0;
        double          error           = 0;
        double          minosError[2]   = { };

        std::unique_ptr<TGraph> upMinScan;
    };

    typedef std::vector<Param> ParamVector;

    int             status      = 0;
    double          chi2        = 0;
    int             ndf         = 0;
    double          prob        = 0;
    double          chi2_ndf    = 0;
    ParamVector     param;

    FitResult() = default;
    FitResult(int status) : status(status) {}
};

////////////////////////////////////////////////////////////////////////////////
FitResult ConstructFitResult( const TFitResult & fitStatus, const ModelCompare::Observable & obs, const FitParamVector & fitParam, int fitIndex,
                              const char * namePrefix, const char * objectiveTitle )
{
    FitResult result;

    result.chi2     = fitStatus.Chi2();
    result.ndf      = fitStatus.Ndf();
    result.prob     = fitStatus.Prob();
    result.chi2_ndf = (result.ndf > 0 ? result.chi2 / result.ndf : 0);

    // get last Fitter object used internally by last Fit call
    TVirtualFitter  *       pFitter      = TVirtualFitter::GetFitter();
    TBackCompFitter *       pBackFitter  = pFitter->InheritsFrom(TBackCompFitter::Class()) ? static_cast<TBackCompFitter *>(pFitter) : nullptr;
    ROOT::Math::Minimizer * pMinimizer   = pBackFitter ? pBackFitter->GetMinimizer() : nullptr;

    for (Int_t i = 0; i < fitParam.size(); ++i)
    {
        if ((fitIndex >= 0) && (i != fitIndex))
            continue;

        FitResult::Param resultParam;

        resultParam.name  = fitParam[i].name;
        resultParam.value = fitStatus.Parameter(i);
        resultParam.error = fitStatus.ParError(i);

        if (fitStatus.HasMinosError(i))
        {
            resultParam.minosError[0] = fitStatus.LowerError(i);
            resultParam.minosError[1] = fitStatus.UpperError(i);
        }

        // create a minimization scan
        if (pMinimizer)
        {
            resultParam.upMinScan.reset( new TGraph( (Int_t)100 ) );
            TGraph * pGraph = resultParam.upMinScan.get();  // alias

            // do +/- 2x error scan (default range)
            double scanMin = resultParam.value - 2 * resultParam.error;
            double scanMax = resultParam.value + 2 * resultParam.error;

            // scan symmetric around 0.0 point
            scanMin = std::min( scanMin, -scanMax );
            scanMax = std::max( scanMax, -scanMin );

            unsigned int nStep = pGraph->GetN();
            bool scanResult = pMinimizer->Scan( i, nStep, pGraph->GetX(), pGraph->GetY(), scanMin, scanMax );
            if (!scanResult || (nStep == 0))
            {
                // failed
                resultParam.upMinScan.reset();
                pGraph = nullptr;
            }
            else
            {
                // success
                pGraph->Set( (Int_t)nStep );  // resize to actual steps

                std::string sName  = std::string(namePrefix) + ((fitIndex < 0) ? "_min_all_"  : "_min_one_") + std::string(obs.name)  + "_" + std::string(resultParam.name);
                std::string sTitle = std::string(obs.title) + ": fit min. for " + std::string(resultParam.name) + ((fitIndex < 0) ? " (all)" : " (one)");

                pGraph->SetName(  sName .c_str() );
                pGraph->SetTitle( sTitle.c_str() );

                pGraph->GetXaxis()->SetTitle( resultParam.name );
                pGraph->GetYaxis()->SetTitle( objectiveTitle   );
            }
        }

        result.param.push_back( std::move(resultParam) );
    }

    return result;
}

////////////////////////////////////////////////////////////////////////////////
class ReweightPDFFunc
{
public:
    typedef std::function<void(const double * par)> TestFunction;

    ReweightPDFFunc(    const FitParamVector & fitParam,
                        const CStringVector & coefNames,
                        const TH1D & source, const ConstTH1DVector & sourceCoefs, const std::vector<double> & sourceEval,
                        double xMin, double xMax,       // used to limit the function to this region
                        double normalization   = 0.0,   // normalizes the reweighted histogram to this value (0.0 skips normalization)
                        double scale           = 1.0,   // F(x) = scale * (content(x))^power
                        double power           = 1.0    // F(x) = scale * (content(x))^power
                    )
    :
        m_nPar(         fitParam.size()     ),
        m_coefNames(    coefNames           ),
        m_source(       source              ),
        m_sourceCoefs(  sourceCoefs         ),
        m_sourceEval(   sourceEval          ),
        m_lastPar(      fitParam.size()     ),
        m_targetParam(  fitParam.size()     ),
        m_xMin(         xMin                ),
        m_xMax(         xMax                ),
        m_norm(         normalization       ),
        m_scale(        scale               ),
        m_power(        power               )
    {
        m_bProfile = source.InheritsFrom(TProfile::Class());

        for (Int_t i = 0; i < m_nPar; ++i)
        {
            m_targetParam[i].name  = fitParam[i].name;
            m_targetParam[i].value = fitParam[i].initValue;
        }
    }

    void EnableThrowOnReject( bool enable )
    {
        m_bThrowOnReject = enable;
    }

    // set optional TestFunction that is called after a new histogram has been constructed
    void SetTestFunction( const TestFunction & func )
    {
        m_testFunc = func;
    }


    Double_t operator()( const Double_t * x, const Double_t * par )
    {
        if (!m_upLastSourceReweight ||
            (memcmp( m_lastPar.data(), par, m_nPar * sizeof(*par)) != 0))
        {
            memcpy( m_lastPar.data(), par, m_nPar * sizeof(*par) );

            std::string sPar;
            for (Int_t i = 0; i < m_nPar; ++i)
            {
                m_targetParam[i].value = par[i];

                if (!sPar.empty()) sPar += ", ";
                sPar += StringFormat( "%g", FMT_F(par[i]) );
            }

            std::vector<double> targetEval;
            CalcEvalVector( m_coefNames, m_targetParam, targetEval );

            m_upLastSourceReweight.reset(
                ReweightEFT::ReweightHist( m_source, m_sourceCoefs, m_sourceEval, targetEval,
                                           "RWFitHist", "RWFitHist" ) );

            double integral         = m_upLastSourceReweight->Integral();
            double integralExt      = ExtendedIntegral( *m_upLastSourceReweight );
            double integralAfter    = integral;
            double integralExtAfter = integralExt;

            if ((m_norm != 0.0) && !m_bProfile)
            {
                double scale = m_norm / integral;
                m_upLastSourceReweight->Scale( scale );

                integralAfter    = m_upLastSourceReweight->Integral();
                integralExtAfter = ExtendedIntegral( *m_upLastSourceReweight );
            }

            //LogMsgInfo( m_norm == 0.0 ? "RWPDFFunc(%hs) integral=%g (%g)" :
            //                            "RWPDFFunc(%hs) integral=%g (%g) -> %g (%g)",
            //            FMT_HS(sPar.c_str()), FMT_F(integral), FMT_F(integralExt),
            //            FMT_F(integralAfter), FMT_F(integralExtAfter) );

            if (m_testFunc)
                m_testFunc( par );
        }

        const TH1D * pHist = m_upLastSourceReweight.get();

        Double_t xVal    = x[0];
        Int_t    bin     = pHist->FindFixBin( xVal );
        Double_t nEff    = GetHistBinEffectiveEntries( *pHist, bin );
        Double_t content = pHist->GetBinContent( bin );

        //LogMsgInfo( "x=%f b=%i n=%f c=%f", FMT_F(xVal), FMT_I(bin), FMT_F(nEff), FMT_F(content) );

        if ((xVal < m_xMin) || (xVal > m_xMax))
        {
            RejectPoint(xVal);
            return 0;
        }

        if ((bin < 1) || (bin > pHist->GetSize() - 2))
        {
            RejectPoint(xVal);
            return 0;
        }

        if (nEff <= 0)
        {
            RejectBin(bin);
            return 0;
        }

        if ((m_power != 1.0) && (m_power != 0.0))
            content = pow( content, m_power );

        if ((m_scale != 1.0) && (m_scale != 0.0))
            content *= m_scale;

        if ((content <= 0.0) && !m_bProfile)
            LogMsgWarning( "RWPDFFunc returning %g", FMT_F(content) );

        return content;
    }


    void Reset()
    {
        m_upLastSourceReweight.reset();
        m_rejectCount = 0;
    }

    size_t RejectCount() const { return m_rejectCount; }

    static double ExtendedIntegral( const TH1D & hist )
    {
        double total = 0.0;

        const Int_t nSize = hist.GetSize();
        for (Int_t bin = 0; bin < nSize; ++bin)
        {
            total += hist.GetBinContent(bin);
        }

        return total;
    }

private:
    void RejectPoint(double x)
    {
        if (m_bThrowOnReject)
            ThrowError( "RWPDFFunc rejecting x=%f", FMT_F(x) );

        LogMsgInfo( "RWPDFFunc rejecting x=%f", FMT_F(x) );
        ++m_rejectCount;
        TF1::RejectPoint();
    }

    void RejectBin(int bin)
    {
        if (m_bThrowOnReject)
            ThrowError( "RWPDFFunc rejecting bin %i", FMT_I(bin) );

        LogMsgInfo( "RWPDFFunc rejecting bin %i", FMT_I(bin) );
        ++m_rejectCount;
        TF1::RejectPoint();
    }

private:
    const size_t                    m_nPar;
    const CStringVector &           m_coefNames;
    const TH1D &                    m_source;
    const ConstTH1DVector &         m_sourceCoefs;
    const std::vector<double> &     m_sourceEval;
    const double                    m_xMin              = -std::numeric_limits<double>::max();
    const double                    m_xMax              =  std::numeric_limits<double>::max();
    const double                    m_norm              = 0.0;
    const double                    m_scale             = 1.0;
    const double                    m_power             = 1.0;
    bool                            m_bProfile          = false;

    std::vector<Double_t>           m_lastPar;
    ReweightEFT::ParamVector        m_targetParam;
    TH1DUniquePtr                   m_upLastSourceReweight;
    size_t                          m_rejectCount       = 0;
    bool                            m_bThrowOnReject    = false;

    std::function<void(const double * par)> m_testFunc  = nullptr;
};

////////////////////////////////////////////////////////////////////////////////
FitResult FitEFTObs( const ModelCompare::Observable & obs, const CStringVector & coefNames,
                     const TH1D & target, const FitParamVector & fitParam,
                     const TH1D & source, const ConstTH1DVector & sourceCoefs, const std::vector<double> & sourceEval,
                     int fitIndex = -1,         // -1 = fit all, otherwise index of fit parameter to fit, keeping all others fixed
                     bool bShapeOnly = false )  // if true normalize event count of the model data to same as the target data
{
    const Double_t  xMin  = target.GetXaxis()->GetXmin();
    const Double_t  xMax  = target.GetXaxis()->GetXmax();
    const Int_t     nPar  = (Int_t)fitParam.size();

    //
    // Construct fit data
    //

    struct FitData
    {
        FitData( const TH1D & data )
        {
            upHist.reset( (TH1D *)data.Clone() );   // clone hist as Fit is not const

            if (upHist->InheritsFrom(TProfile::Class()))
            {
                TProfile * pProfile = static_cast<TProfile *>(upHist.get());
                pProfile->SetErrorOption("");
                upGraph.reset( GraphFromProfile( *pProfile ) );
            }
        }

        TFitResultPtr Fit( TF1 * f1, Option_t * option = "" )
        {
            if (upGraph)
                return upGraph->Fit( f1, (std::string(option) + " EX0").c_str() );
            else
                return upHist ->Fit( f1, option );
        }

        Double_t Chisquare( TF1 * f1 )
        {
            if (upGraph)
                return upGraph->Chisquare( f1 );
            else
                return upHist ->Chisquare( f1 );
        }

        bool IsHist()    const { return !upGraph;  }
        bool IsProfile() const { return !!upGraph; }

        std::unique_ptr<TH1D>       upHist;
        std::unique_ptr<TGraph>     upGraph;
    };

    FitData fitData( target );

    bool bLogLike = fitData.IsHist();

    //
    // setup fit model function
    //

    const double normalization = !(bShapeOnly && bLogLike) ? 0.0 : target.Integral(); // ReweightPDFFunc::ExtendedIntegral(target);

    ReweightPDFFunc modelFunc(  fitParam, coefNames, source, sourceCoefs, sourceEval,
                                xMin, xMax, normalization );

    TF1 fitFunc( "EFT", &modelFunc, xMin, xMax, nPar );

    for (Int_t i = 0; i < nPar; ++i)
    {
        fitFunc.SetParName(   i, fitParam[i].name );
        fitFunc.SetParameter( i, fitParam[i].initValue );

        if ((fitIndex >= 0) && (i != fitIndex))
        {
            fitFunc.FixParameter( i, fitParam[i].initValue );
            fitFunc.SetParError(  i, 0.0 );
        }
        else
        {
            fitFunc.SetParLimits( i, fitParam[i].minValue, fitParam[i].maxValue );
            fitFunc.SetParError(  i, (fitParam[i].maxValue - fitParam[i].minValue) / 1E6 ); // set initial step (helps fit converge)
        }
    }

    //
    // setup fit options
    //

    // do not use M - it produces too many warnings of:
    // "FUNCTION VALUE DOES NOT SEEM TO DEPEND ON ANY OF THE 1 VARIABLE PARAMETERS. VERIFY THAT STEP SIZES ARE BIG ENOUGH AND CHECK FCN LOGIC."

    std::string fitOption1 = "N S Q";
    std::string fitOption2 = "N S E";  // E=hesse error (invert error matrix) and run minos errors

    // use log-likelihood for TH1D but not TProfile
    // Note: do not use "WL" as it does something weird

    if (bLogLike)
    {
        fitOption1 += " L";
        fitOption2 += " L";
    }

    //
    // Setup minimizer options
    //

    RestoreDefaultMinimizerOptions();

    // increase error definition to 2 sigma
    // error definition: the objective function delta that determines parameter error results
    // default is 1 sigma - likelihood: 0.5  chi^2: 1.0
    // both have quadratic dependence, so 2 sigma is 4 times the 1 sigma error definition
    ROOT::Math::MinimizerOptions::SetDefaultErrorDef( bLogLike ? 2 : 4 );

    //
    // first fit with limits
    //

    // clear last fitter to ensure no previous fit results will affect this fit
    // should not be necessary, but paranoia avoids issues
    ClearLastFitter();

    TFitResultPtr fitStatus;

    fitStatus = fitData.Fit( &fitFunc, fitOption1.c_str() );
    if ((int)fitStatus != 0)
    {
        fitStatus->Print();
        if (((int)fitStatus < 0) || ((int)fitStatus % 1000 != 0))  // ignore improve (M) errors
            return FitResult( (int)fitStatus );
    }

    //LogMsgInfo( "-*-*-*-*-*-" );

    //
    // second fit with no limits
    //

    // remove parameter limits on non-fixed parameters
    for (Int_t i = 0; i < nPar; ++i)
    {
        if ((fitIndex < 0) || (i == fitIndex))
            fitFunc.ReleaseParameter(i);
    }

    modelFunc.Reset();

    // re-run fit using result from previous fit as initial values
    // including parameter errors that are used to setup fit steps
    fitStatus = fitData.Fit( &fitFunc, fitOption2.c_str() );
    if ((int)fitStatus != 0)
    {
        fitStatus->Print();
        if (((int)fitStatus < 0) || ((int)fitStatus % 1000 != 0))   // ignore improve (M) errors
            return FitResult( (int)fitStatus );
    }

    /*
    if (pFitHist->InheritsFrom(TProfile::Class()))
    {
        // TEST: What happens if I double the errors?
        // Result: - doubles the error on cWWWW, cW, cB
        //         - reduces chi2 error by 1/4

        TProfile * pProf  = static_cast<TProfile *>(pFitHist.get());
        Double_t * pSumw2 = pProf->GetSumw2()->GetArray();
        if (pSumw2)
        {
            //LogMsgHistDump(*pProf);

            const Int_t nSize = pProf->GetSize();
            for (Int_t bin = 0; bin < nSize; ++bin)
            {
                if (bin == 55)
                    bin = bin;

                Double_t binContent = pProf->GetBinContent(bin);
                Double_t binEntries = pProf->GetBinEntries(bin);
                pSumw2[bin] = 4 * pSumw2[bin] - 3 * binEntries * binContent * binContent;
            }

            //LogMsgHistDump(*pProf);

            fitStatus = pFitHist->Fit( &fitFunc, fitOption2.c_str() );
            if ((int)fitStatus != 0)
            {
                if (((int)fitStatus < 0) || ((int)fitStatus % 1000 != 0))   // ignore improve (M) errors
                    return FitResult( (int)fitStatus );
            }
        }
    }
    */

    fitStatus->Print();

    LogMsgInfo( "\nReject Count: %u", FMT_U(modelFunc.RejectCount()) );

    // fill in result

    FitResult result = ConstructFitResult( *fitStatus, obs, fitParam, fitIndex, "binned", bLogLike ? "Log likelihood" : "#chi^{2}" );

    // cross-check chi2
    {
        Double_t chi2 = fitData.Chisquare( &fitFunc );

        LogMsgInfo( "Cross-check chi2: %g", FMT_F(chi2) );
    }

    LogMsgHistBinCounts( target );

    LogMsgInfo( "" );  // empty line

    return result;
}

////////////////////////////////////////////////////////////////////////////////
FitResult UnBinFitEFTObs( const ModelCompare::Observable & obs, const CStringVector & coefNames, const FitParamVector & fitParam,
                          const TNtupleD & target, double targetScale,
                          const TH1D & source, const ConstTH1DVector & sourceCoefs, const std::vector<double> & sourceEval,
                          int fitIndex = -1 )   // -1 = fit all, otherwise index of fit parameter to fit, keeping all others fixed
{
    /*
    static int count = 0;

    if (count++ >= 8)
        return FitResult();
    */

    const Double_t  xMin  = source.GetXaxis()->GetXmin();
    const Double_t  xMax  = source.GetXaxis()->GetXmax();
    const Int_t     nPar  = (Int_t)fitParam.size();

    //
    // setup fit model function
    //

    TF1 * pFitFunc = nullptr;

    /*
    std::unique_ptr<ROOT::Fit::UnBinData> upTestData = nullptr;

    auto testFunc = [&](const double * par) -> void
    {
        auto func = ROOT::Math::WrappedMultiTF1(*pFitFunc, 1);

        unsigned int nPoints(0);
        double logL = ROOT::Fit::FitUtil::EvaluateLogL( func, *upTestData, par, 0, false, nPoints );


        const ROOT::Fit::UnBinData & test = *upTestData;

        size_t skipped = 0;
        double myLogL = 0;
        unsigned int n = test.Size();
        for (unsigned int i = 0; i < n; ++i)
        {
            const double * x = test.Coords(i);
            double f = func( x, par );
            if (f <= 0.0)
            {
                ++skipped;
                continue;
            }
            myLogL += ROOT::Math::Util::EvalLog(f);
        }
        myLogL = -myLogL;

        if (logL == myLogL)
            LogMsgInfo( "LogL = %.15e (skipped=%u)\n", FMT_F(logL), FMT_U(skipped) );
        else
            LogMsgInfo( "LogL = %.15e myLogL = %.15e (skipped=%u)\n", FMT_F(logL), FMT_F(myLogL), FMT_U(skipped) );

    };
    */

    // Create model function
    //
    // For negative log-likelihood, the scale and power parameters passed to
    // the ReweightPDFFunc constructor affect the objective as so:
    //   FCN = power * -log(content) - log(scale)
    //
    // The PDF must be normalized to some constant for the fit to function, unless
    // one performs an extended log-likelhood fit (fitOption.Like = 1).
    // An extended log-likelihood fit adds a Poisson term, that adds the integral
    // of the PDF (in the range) to the result. See FitUtil.cxx, EvaluateLogL().
    // This adds whatever we normalize to, addng 1 to the result.

    ReweightPDFFunc modelFunc(  fitParam, coefNames, source, sourceCoefs, sourceEval, xMin, xMax,
                                1.0,            // normalize PDF to 1.0 (required for unbinned fits)
                                1.0,            // default scale
                                //1.0           // default power
                                targetScale     // multiply log-likelood by targetScale
                              );

    //modelFunc.SetTestFunction( testFunc );

    pFitFunc = new TF1( "EFT_unbinned", &modelFunc, xMin, xMax, nPar );
    TF1 & fitFunc = *pFitFunc;

    for (Int_t i = 0; i < nPar; ++i)
    {
        fitFunc.SetParName(   i, fitParam[i].name );
        fitFunc.SetParameter( i, fitParam[i].initValue );

        if ((fitIndex >= 0) && (i != fitIndex))
        {
            fitFunc.FixParameter( i, fitParam[i].initValue );
            fitFunc.SetParError(  i, 0.0 );
        }
        else
        {
            fitFunc.SetParLimits( i, fitParam[i].minValue, fitParam[i].maxValue );
            //fitFunc.SetParError(  i, (fitParam[i].maxValue - fitParam[i].minValue) / 1E6 ); // set initial step (helps fit converge)
        }
    }

    //
    // Construct fit data
    //

    std::vector<double> data;
    {
        std::unique_ptr<TNtupleD> upRead( (TNtupleD *)target.Clone() );

        Long64_t nEntries = upRead->GetEntries();

        data.reserve( (size_t)nEntries );

        //upTestData.reset( new ROOT::Fit::UnBinData( (unsigned int)nEntries ) );

        for (Long64_t entry = 0; entry < nEntries; ++entry)
        {
          //upRead->LoadTree(entry);    // non-const
            upRead->GetEntry(entry);    // non-const

            Double_t value = *upRead->GetArgs();

            if ((value < xMin) || !(value < xMax))
                continue;

            // omit points rejected by our fit function
            TF1::RejectPoint(false);
            fitFunc( &value ); // evaluate using stored function parameters
            if (TF1::RejectedPoint())
                continue;

            data.push_back( value );

            //upTestData->Add( value );
        }
    }

    //
    // setup fit options
    //

    Foption_t fitOption;
    {
        //fitOption.More          = 1;
        fitOption.StoreResult   = 1;
        fitOption.Nograph       = 1;
        fitOption.Verbose       = 1;
      //fitOption.Errors        = 1;
      //fitOption.Like          = 1;  // 1 = extended (adds Poisson term). Requires range to be set in the unbinned data object.
    }

    //
    // Setup minimizer options
    //

    RestoreDefaultMinimizerOptions();

    ROOT::Math::MinimizerOptions minOption;
    {
        // for likelihood: 1 sigma = 0.5*1^2 = 0.5, 2 sigma = 0.5*2^2 = 2
        minOption.SetErrorDef( 2.0 );
        //minOption.SetErrorDef( 2.0 / targetScale );

        // converge when EDM < 0.001 * Tolerance
        // default: Tolerance = 0.1
        minOption.SetTolerance( 1000 );

        //minOption.SetPrecision( std::numeric_limits<double>::min() );
        //minOption.SetPrecision( std::numeric_limits<double>::epsilon() );
        //minOption.SetPrecision( 1E-100 );

        minOption.SetStrategy(2);
    }

    //
    // first fit with limits
    //

    // enable throw on reject to check that our PDF always has a value for every unbinned data point
    // unsupported data points should have been removed above as the unbinned data was constructed
    modelFunc.EnableThrowOnReject(true);

    // clear last fitter to ensure no previous fit results will affect this fit
    // should not be necessary, but paranoia avoids issues
    ClearLastFitter();

    ROOT::Fit::DataRange range( xMin, xMax );

    ROOT::Fit::UnBinData *  pFitData = nullptr;
    TFitResultPtr           fitStatus;

    pFitData  = new ROOT::Fit::UnBinData( (unsigned int)data.size(), data.data(), range ); // copies data
    
    fitStatus = ROOT::Fit::UnBinFit( pFitData, &fitFunc, fitOption, minOption );    // TVirtualFitter::GetFitter() takes ownership of pFitData
    if ((int)fitStatus != 0)
    {
        fitStatus->Print();
        if (((int)fitStatus < 0) || ((int)fitStatus % 1000 != 0))  // ignore improve (M) errors
            return FitResult( (int)fitStatus );
    }

    /*
    //
    // second fit with no limits
    //

    // remove parameter limits on non-fixed parameters
    for (Int_t i = 0; i < nPar; ++i)
    {
        if ((fitIndex < 0) || (i == fitIndex))
            fitFunc.ReleaseParameter(i);
    }

    modelFunc.Reset();

    // re-run fit using result from previous fit as initial values
    // including parameter errors that are used to setup fit steps

    pFitData  = new ROOT::Fit::UnBinData( (unsigned int)data.size(), data.data(), range ); // copies data
    
    fitStatus = ROOT::Fit::UnBinFit( pFitData, &fitFunc, fitOption, minOption );    // TVirtualFitter::GetFitter() takes ownership of pFitData
    if ((int)fitStatus != 0)
    {
        fitStatus->Print();
        if (((int)fitStatus < 0) || ((int)fitStatus % 1000 != 0))  // ignore improve (M) errors
            return FitResult( (int)fitStatus );
    }
    */

//    fitStatus->Print();

    LogMsgInfo( "\nReject Count: %u", FMT_U(modelFunc.RejectCount()) );

    // fill in result

    FitResult result = ConstructFitResult( *fitStatus, obs, fitParam, fitIndex, "Unbinned", "Log likelihood" );

    LogMsgInfo( "" );  // empty line

    return result;
}

////////////////////////////////////////////////////////////////////////////////
static TH1D * CreateReweightFitResult( const CStringVector & coefNames,
                                       const FitResult & fitResult, const FitParamVector & fitParam,
                                       const TH1D & source, const ConstTH1DVector & sourceCoefs, const std::vector<double> & sourceEval )
{
    // construct target parameter values
    ReweightEFT::ParamVector targetParam;
    {
        for (const auto & p : fitParam)
        {
            ReweightEFT::Parameter t;
            t.name  = p.name;
            t.value = p.initValue;
            targetParam.push_back(t);
        }

        for (const auto & p : fitResult.param)
        {
            auto itr = std::find_if( targetParam.begin(), targetParam.end(), [&p](const ReweightEFT::Parameter & x) { return strcmp(x.name, p.name) == 0; } );
            if (itr == targetParam.end())
                ThrowError( "Cannot find fit result parameter " + std::string(p.name) );
            itr->value = p.value;
        }

        //for (const auto & t : targetParam)
        //{
        //    LogMsgInfo( "> %hs = %g", FMT_HS(t.name), FMT_F(t.value) );
        //}
    }

    // construct target evaluation
    std::vector<double> targetEval;
    CalcEvalVector( coefNames, targetParam, targetEval );

    TH1D * pHist = ReweightEFT::ReweightHist( source, sourceCoefs, sourceEval, targetEval, "RWFitHist", "RWFitHist" );

    return pHist;
}

////////////////////////////////////////////////////////////////////////////////
static void CrossCheckFitResult( const CStringVector & coefNames,
                                 const FitResult & fitResult, const FitParamVector & fitParam,
                                 const TH1D & target,
                                 const TH1D & source, const ConstTH1DVector & sourceCoefs, const std::vector<double> & sourceEval )
{
    if (fitResult.status != 0)
        return;

    LogMsgInfo( "" );

    TH1DUniquePtr upFitHist
    {
        CreateReweightFitResult( coefNames, fitResult, fitParam,
                                 source, sourceCoefs, sourceEval )
    };

    LogMsgHistEffectiveEntries( target );
    LogMsgHistEffectiveEntries( *upFitHist );

    LogMsgHistBinCounts( target );
    LogMsgHistBinCounts( *upFitHist );

    LogMsgInfo( "" );

    // perform chi2 test
    // note: this test includes the errors from both histograms and thus will not match
    // the chi2 result from the fit, which only uses the target errors.
    // This chi2 value should be lower than the fit chi2.
    
    RootUtil::Chi2Result chi2;
    chi2.Chi2Test( target, *upFitHist );  // supports both TH1D and TProfile

    LogMsgInfo( chi2.Label().c_str() );
}

////////////////////////////////////////////////////////////////////////////////
void LoadTupleData( const ModelCompare::ObservableVector & observables, const ModelCompare::ModelFile & model,
                    TupleVector & modelData,
                    const char * cacheFileName )
{
    modelData.clear();

    bool        bLoadEvents = false;
    TupleVector loadData;

    for (const ModelCompare::Observable & obs : observables)
    {
        if (obs.nDim != 1)
        {
            loadData .push_back(nullptr);
            modelData.push_back(nullptr);
            continue;
        }

        // make a n-tuple
        TNtupleD * pTuple = nullptr;
        {
            std::string sName  = "data_" + obs.BuildHistName(  model.modelName  );
            std::string sTitle =           obs.BuildHistTitle( model.modelTitle );

            pTuple = new TNtupleD( sName.c_str(), sTitle.c_str(), obs.name );
            pTuple->SetDirectory( nullptr );  // decouple from any directory
        }

        if (LoadCacheTuple( cacheFileName, pTuple ))
        {
            LogMsgInfo( "Loaded %hs from cache", FMT_HS(pTuple->GetName()) );
            loadData.push_back( nullptr );  // skip this histogram
        }
        else
        {
            loadData.push_back( pTuple );
            bLoadEvents = true;
        }

        modelData.push_back(pTuple);
    }

    // fill the n-tuple

    auto FillFunc = [&](const HepMC::GenVertex & signal)
    {
        size_t obsIndex = 0;
        for (const ModelCompare::Observable & obs : observables)
        {
            TNtupleD * pTuple = loadData[obsIndex++];
            if (pTuple)
            {
                double value(0);
                obs.getFunction( signal, &value, 1 );
                pTuple->Fill(&value);
            }
        }
    };

    if (bLoadEvents)
    {
        LoadEvents( model.fileName, FillFunc );

        SaveTuples( cacheFileName, ToConstTupleVector(modelData) );
    }

    // debug
    //SaveTuples( "optbin/Cross_Unbinned.root", ToConstTupleVector(modelData) );
}

////////////////////////////////////////////////////////////////////////////////
void LoadTupleData( const ModelCompare::ObservableVector & observables, const ModelCompare::ModelFileVector & models,
                    std::vector<TupleVector> & allData,
                    const char * cacheFileName )
{
    allData.clear();

    for (const ModelCompare::ModelFile & model : models)
    {
        TupleVector modelData;
        LoadTupleData( observables, model, modelData, cacheFileName );
        allData.push_back( modelData );
    }
}

////////////////////////////////////////////////////////////////////////////////
static void WriteGraph( TFile & file, const TGraph & graph )
{
    TDirectory * oldDir = gDirectory;
    if (gDirectory != &file)
        file.cd();

    graph.Write( nullptr, TObject::kOverwrite );

    if (gDirectory != oldDir)
        oldDir->cd();
}

////////////////////////////////////////////////////////////////////////////////
static void WriteLog( FILE * fpLog, const char * format, ... )
{
    std::va_list args;

    va_start(args, format);
    LogMsg( nullptr, format, args, stdout );
    va_end(args);

    va_start(args, format);
    LogMsg( nullptr, format, args, fpLog );
    va_end(args);
}

////////////////////////////////////////////////////////////////////////////////
static bool WriteFitResult( FILE * fpLog, const FitResult & result, TFile * pFigureFile = nullptr )
{
    if (result.status != 0)
    {
        WriteLog( fpLog, "!!! Fit Failed: status = %i", FMT_I(result.status) );
        return false;
    }

    WriteLog( fpLog, "----------------------------------------------------------------------------------------------------" );

    WriteLog( fpLog, "chi2 / ndf = %g / %g = %g (%.2g)", FMT_F(result.chi2), FMT_F(result.ndf), FMT_F(result.chi2_ndf), FMT_F(result.chi2_ndf) );
    WriteLog( fpLog, "prob = %g\n", FMT_F(result.prob) );

    for (const FitResult::Param & p : result.param)
    {
        WriteLog( fpLog, "%8hs = %E ± %E [%.2E, %.2E] | %8.2g ± %-6.2g [%.2g, %.2g]  E-6",
                  FMT_HS(p.name),
                  FMT_F(p.value),     FMT_F(p.error),     FMT_F(p.minosError[0]),     FMT_F(p.minosError[1]),
                  FMT_F(p.value*1E6), FMT_F(p.error*1E6), FMT_F(p.minosError[0]*1E6), FMT_F(p.minosError[1]*1E6)
                );

        if (pFigureFile && p.upMinScan)
            WriteGraph( *pFigureFile, *p.upMinScan );
    }

    WriteLog( fpLog, "----------------------------------------------------------------------------------------------------" );

    return true;
}

////////////////////////////////////////////////////////////////////////////////
void FitEFT( const char * outputFileName,
             const ModelCompare::ObservableVector & observables,
             const CStringVector & coefNames,
             const ModelCompare::ModelFile & targetFile, const FitParamVector & fitParam,
             const ModelCompare::ModelFile & sourceFile, const ReweightEFT::ParamVector & sourceParam,
             double luminosity,
             const char * cacheFileName )
{
    // disable automatic histogram addition to current directory
    TH1::AddDirectory(kFALSE);
    // enable automatic sumw2 for every histogram
    TH1::SetDefaultSumw2(kTRUE);

    // disable reusing the same minuit object
    // this fixes the observed behavior/bug where a fit affected subsequent fits
    TMinuitMinimizer::UseStaticMinuit(false);

    // ------

    std::string outFileName = std::string(outputFileName) + ".root";
    std::string logFileName = std::string(outputFileName) + ".txt";

    std::string cacheHistFileName(cacheFileName);
    std::string cacheTreeFileName(cacheFileName);

    if (!cacheHistFileName.empty()) cacheHistFileName += "_hist.root";
    if (!cacheTreeFileName.empty()) cacheTreeFileName += "_tree.root";

    LogMsgInfo( "Root output file: %hs", FMT_HS(outFileName.c_str()) );
    std::unique_ptr<TFile> upOutputFile( new TFile( outFileName.c_str(), "RECREATE" ) );
    if (upOutputFile->IsZombie() || !upOutputFile->IsOpen())    // IsZombie is true if constructor failed
    {
        LogMsgError( "Failed to create output file (%hs).", FMT_HS(outFileName.c_str()) );
        ThrowError( std::invalid_argument( outFileName ) );
    }

    LogMsgInfo( "Log output file: %hs", FMT_HS(logFileName.c_str()) );
    FILE * fpLog = fopen( logFileName.c_str(), "wt" );

    // calculate luminosity scales

    double targetScale = luminosity * targetFile.crossSection * 1000 / targetFile.crossSectionEvents;

    // load input trees

    TupleVector targetTrees;    // targetTrees[observable]

    LoadTupleData( observables, targetFile, targetTrees, cacheTreeFileName.c_str() );

    // load input histograms

    TH1DVector              targetHists;        // targetHists[observable]
    TH1DVector              sourceHists;        // sourceHists[observable]
    std::vector<TH1DVector> sourceCoefs;        // sourceCoefs[observable][coefficient]
    std::vector<double>     sourceEval;         // sourceEval[coefficient]
    TH1DVector              rawtargetHists;     // rawtargetHists[observable]
    TH1DVector              rawsourceHists;     // rawsourceHists[observable]

    ReweightEFT::LoadReweightFiles( observables, coefNames, targetFile, sourceFile, sourceParam,    // inputs
                                    targetHists, sourceHists, sourceCoefs, sourceEval,              // outputs
                                    rawtargetHists, rawsourceHists,
                                    luminosity, cacheHistFileName.c_str() );                        // optionals


    // write input hists

    WriteHists( upOutputFile.get(), targetHists );   // output file takes ownership of histograms
    WriteHists( upOutputFile.get(), sourceHists );   // output file takes ownership of histograms

    // uncomment to save coefficent hists to file
    //  for (const auto & coefData : sourceCoefs)
    //      WriteHists( upOutputFile.get(), coefData );  // output file takes ownership of histograms

    LogMsgInfo("");
    LogMsgHistBinCounts( ToConstTH1DVector(targetHists) );
    LogMsgInfo("");
    LogMsgHistBinCounts( ToConstTH1DVector(sourceHists) );
    LogMsgInfo("");
    LogMsgHistBinCounts( ToConstTH1DVector(targetHists), ToConstTH1DVector(sourceHists) );
    LogMsgInfo("");

    /*
    // replace targetHists errors with sqrt(N) errors - Asimov samples
    // The current errors on targetHists, are errors on the expectation value.
    // To convert this to pseudo-data, we keep the expectation value, but replace
    // the errors with sqrt(binContent).

    LogMsgInfo("------------ Before Asimov ------------");
    targetHists[0]->Print("all");
    LogMsgHistStats(*targetHists[0]);

    for (TH1D * pHist : targetHists)
    {
        pHist->Sumw2(kFALSE);
        pHist->Sumw2(kTRUE);
        pHist->ResetStats();
    }

    LogMsgInfo("------------ After Asimov ------------");
    LogMsgHistStats(*targetHists[0]);
    targetHists[0]->Print("all");
    */

    // fit each source hist to the target, varying only one parameter

    typedef std::tuple<const char *, bool, int>     FitIdent;           // obsName, binned, fitIndex
    typedef std::map< FitIdent, FitResult >         FitResultMap;

    FitResultMap results;

    size_t fitFail = 0;

    for (size_t obsIndex = 0; obsIndex < observables.size(); ++obsIndex)
    {
        const ModelCompare::Observable &    obs         = observables[obsIndex];
        const TNtupleD *                    pTargetTree = targetTrees[obsIndex];
        const TH1D *                        pTargetHist = targetHists[obsIndex];
        const TH1D &                        source      = *sourceHists[obsIndex];
        const ConstTH1DVector               srcCoefs    = ToConstTH1DVector(sourceCoefs[obsIndex]);

        ModelCompare::GoodBadHists goodBad;

        WriteLog( fpLog, "\n******************** %hs ********************", FMT_HS(obs.name) );

        if (pTargetHist->InheritsFrom(TProfile::Class()))
        {
            goodBad = ModelCompare::HistSplitGoodBadBins( pTargetHist,        rawtargetHists[obsIndex] );
            goodBad = ModelCompare::HistSplitGoodBadBins( goodBad.good.get(), rawsourceHists[obsIndex] );
            pTargetHist = goodBad.good.get();
        }

        for (int fitIndex = 0; fitIndex < fitParam.size(); ++fitIndex)
        {
            if (pTargetTree)
            {
                WriteLog( fpLog, "\n------------------------------------------------------------" );
                WriteLog( fpLog, "Fitting unbinned %hs to parameter %hs",
                                    FMT_HS(obs.name), FMT_HS(fitParam[fitIndex].name) );
                LogMsgInfo( "------------------------------------------------------------" );

                FitResult fitResult = UnBinFitEFTObs( obs, coefNames, fitParam, *pTargetTree, targetScale,
                                                      source, srcCoefs, sourceEval, fitIndex );

                if (!WriteFitResult( fpLog, fitResult, upOutputFile.get() ))
                    ++fitFail;

                //exit(1);

                results[ FitIdent(obs.name,false,fitIndex) ] = std::move(fitResult);
            }

            if (pTargetHist)
            {
                WriteLog( fpLog, "\n------------------------------------------------------------" );
                WriteLog( fpLog, "Fitting binned %hs to parameter %hs",
                                    FMT_HS(obs.name), FMT_HS(fitParam[fitIndex].name) );
                LogMsgInfo( "------------------------------------------------------------" );

                FitResult fitResult = FitEFTObs( obs, coefNames, *pTargetHist, fitParam,
                                                 source, srcCoefs, sourceEval, fitIndex );

                if (!WriteFitResult( fpLog, fitResult, upOutputFile.get() ))
                    ++fitFail;

                // cross-check
                CrossCheckFitResult( coefNames, fitResult, fitParam, *pTargetHist,
                                     source, srcCoefs, sourceEval );

                results[ FitIdent(obs.name,true,fitIndex) ] = std::move(fitResult);
            }
        }

        if (pTargetTree)
        {
            WriteLog( fpLog, "\n------------------------------------------------------------" );
            WriteLog( fpLog, "Fitting unbinned %hs to all parameters",
                                FMT_HS(obs.name) );
            LogMsgInfo( "------------------------------------------------------------" );

            FitResult fitResult = UnBinFitEFTObs( obs, coefNames, fitParam, *pTargetTree, targetScale,
                                                  source, srcCoefs, sourceEval, -1 );

            if (!WriteFitResult( fpLog, fitResult, upOutputFile.get() ))
                ++fitFail;

            results[ FitIdent(obs.name,false,-1) ] = std::move(fitResult);
        }

        if (pTargetHist)
        {
            WriteLog( fpLog, "\n------------------------------------------------------------" );
            WriteLog( fpLog, "Fitting binned %hs to all parameters",
                                FMT_HS(obs.name) );
            LogMsgInfo( "------------------------------------------------------------") ;

            FitResult fitResult = FitEFTObs( obs, coefNames, *pTargetHist, fitParam,
                                             source, srcCoefs, sourceEval, -1 );

            if (!WriteFitResult( fpLog, fitResult, upOutputFile.get() ))
                ++fitFail;

            // cross-check
            CrossCheckFitResult( coefNames, fitResult, fitParam, *pTargetHist,
                                 source, srcCoefs, sourceEval );

            results[ FitIdent(obs.name,true,-1) ] = std::move(fitResult);
        }
    }

    // TODO: make figures from fit result

    WriteLog( fpLog, "\nSUMMARY" );

    // Write summary
    for (int fitIndex = 0; fitIndex < fitParam.size(); ++fitIndex)
    {
        const FitParam & fpar = fitParam[fitIndex];

        WriteLog( fpLog, "\n---------- %hs ----------", FMT_HS(fpar.name) );

        for (size_t obsIndex = 0; obsIndex < observables.size(); ++obsIndex)
        {
            const ModelCompare::Observable & obs = observables[obsIndex];

            bool bProfile = sourceHists[obsIndex]->InheritsFrom(TProfile::Class());

            for (size_t binnedFlag = 0; binnedFlag <= 1; ++binnedFlag)
            {
                bool bBinned = !!binnedFlag;

                if (bProfile && !bBinned)
                    continue;

                // find fit of one parameter

                const FitResult *        pFitOne = nullptr;
                const FitResult::Param * pParOne = nullptr;
                {
                    auto itrRes = results.find( FitIdent(obs.name,bBinned,fitIndex) );
                    if (itrRes != results.end())
                    {
                        pFitOne = &itrRes->second;

                        auto itrPar = std::find_if( pFitOne->param.cbegin(), pFitOne->param.cend(), [&](const FitResult::Param & p) { return strcmp(p.name, fpar.name) == 0; } );
                        if (itrPar != pFitOne->param.cend())
                            pParOne = &*itrPar;
                    }
                }

                // find fit of all parameters

                const FitResult *        pFitAll = nullptr;
                const FitResult::Param * pParAll = nullptr;
                {
                    auto itrRes = results.find( FitIdent(obs.name,bBinned,-1) );
                    if (itrRes != results.end())
                    {
                        pFitAll = &itrRes->second;

                        auto itrPar = std::find_if( pFitAll->param.cbegin(), pFitAll->param.cend(), [&](const FitResult::Param & p) { return strcmp(p.name, fpar.name) == 0; } );
                        if (itrPar != pFitAll->param.cend())
                            pParAll = &*itrPar;
                    }
                }

                auto GetResult = [](const FitResult::Param * pPar, const FitResult * pRes, double scale = 1E6) -> std::string
                {
                    auto GetMinos = [](const FitResult::Param & par, double scale) -> std::string
                    {
                        double err = par.error         * scale;
                        double lwr = par.minosError[0] * scale;
                        double upr = par.minosError[1] * scale;

                        if ((lwr == 0.0) && (upr == 0.0)) return std::string();

                        std::string sErr = StringFormat( "%.2g", FMT_F(err)  );
                        std::string sLwr = StringFormat( "%.2g", FMT_F(-lwr) );
                        std::string sUpr = StringFormat( "%.2g", FMT_F(upr)  );

                        if ((sLwr == sErr) && (sUpr == sErr)) return std::string();

                        return StringFormat( "(%+.2g,%+.2g)", FMT_F(lwr), FMT_F(upr) );
                    };

                    if (!pPar)
                        return StringFormat( "%8hs | %-28hs | %6hs", "", "Failed", "" );

                    std::string sMinos  = GetMinos( *pPar, scale );
                    std::string sChi2   = (!pRes || pRes->chi2_ndf < 0.0) ? "" : StringFormat( "%6.2g", FMT_F(pRes->chi2_ndf) );

                    return StringFormat( "%8.2g | %8.2g %19hs | %6hs",
                                         FMT_F(pPar->value*scale), FMT_F(pPar->error*scale),
                                         FMT_HS(sMinos.c_str()), FMT_HS(sChi2.c_str()) );
                };

                std::string sParOne = GetResult( pParOne, pFitOne );
                std::string sParAll = GetResult( pParAll, pFitAll );

                WriteLog( fpLog, "%C %12hs: %48s || %48s",
                          FMT_HC(bBinned ? 'B' : 'U'), FMT_HS(obs.name),
                          FMT_HS(sParOne.c_str()), FMT_HS(sParAll.c_str()) );
            }
        }
    }

    if (fitFail)
        WriteLog( fpLog, "\n!!! Error: %u fit(s) failed !!!", FMT_U(fitFail) );

    WriteLog( fpLog, "" );

    fclose(fpLog);
}

////////////////////////////////////////////////////////////////////////////////

}  // namespace FitEFT
