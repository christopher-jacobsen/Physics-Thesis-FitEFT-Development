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

using namespace RootUtil;

////////////////////////////////////////////////////////////////////////////////
namespace FitEFT
{

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
class ReweightPDFFunc
{
public:
    ReweightPDFFunc( const CStringVector & coefNames,
                     const FitParamVector & fitParam,
                     const TH1D & source, const ConstTH1DVector & sourceCoefs, const std::vector<double> & sourceEval )
    :
        m_nPar(        fitParam.size()  ),
        m_coefNames(   coefNames        ),
        m_source(      source           ),
        m_sourceCoefs( sourceCoefs      ),
        m_sourceEval(  sourceEval       ),
        m_lastPar(     fitParam.size()  ),
        m_targetParam( fitParam.size()  )
    {
        for (Int_t i = 0; i < m_nPar; ++i)
        {
            m_targetParam[i].name  = fitParam[i].name;
            m_targetParam[i].value = fitParam[i].initValue;
        }
    }

    Double_t operator()( const Double_t * x, const Double_t * par )
    {
        if (!m_upLastSourceReweight ||
            (memcmp( m_lastPar.data(), par, m_nPar * sizeof(*par)) != 0))
        {
            for (Int_t i = 0; i < m_nPar; ++i)
                m_targetParam[i].value = par[i];

            std::vector<double> targetEval;
            CalcEvalVector( m_coefNames, m_targetParam, targetEval );

            m_upLastSourceReweight.reset(
                ReweightEFT::ReweightHist( m_source, m_sourceCoefs, m_sourceEval, targetEval,
                                           "RWFitHist", "RWFitHist" ) );

            memcpy( m_lastPar.data(), par, m_nPar * sizeof(*par) );
        }

        const TH1D * pHist = m_upLastSourceReweight.get();

        Double_t xVal    = x[0];
        Int_t    bin     = pHist->FindFixBin( xVal );
        Double_t nEff    = GetHistBinEffectiveEntries( *pHist, bin );
        Double_t content = pHist->GetBinContent( bin );

        //LogMsgInfo( "x=%f b=%i n=%f c=%f", FMT_F(xVal), FMT_I(bin), FMT_F(nEff), FMT_F(content) );

        // reduce NDF
        if (nEff <= 0)
        {
            // LogMsgInfo( ">>>> Rejecting bin %i <<<<", bin );
            ++m_rejectCount;
            TF1::RejectPoint();
            return 0;
        }

        return content;
    }

    void Reset()
    {
        m_upLastSourceReweight.reset();
        m_rejectCount = 0;
    }

    size_t RejectCount() const { return m_rejectCount; }

private:
    const size_t                    m_nPar;
    const CStringVector &           m_coefNames;
    const TH1D &                    m_source;
    const ConstTH1DVector &         m_sourceCoefs;
    const std::vector<double> &     m_sourceEval;

    std::vector<Double_t>           m_lastPar;
    ReweightEFT::ParamVector        m_targetParam;
    ConstTH1DUniquePtr              m_upLastSourceReweight;
    size_t                          m_rejectCount               = 0;
};

////////////////////////////////////////////////////////////////////////////////
FitResult FitEFTObs( const ModelCompare::Observable & obs, const CStringVector & coefNames,
                     const TH1D & target, const FitParamVector & fitParam,
                     const TH1D & source, const ConstTH1DVector & sourceCoefs, const std::vector<double> & sourceEval,
                     int fitIndex = -1 )   // -1 = fit all, otherwise index of fit parameter to fit, keeping all others fixed
{
    const Double_t  xmin  = target.GetXaxis()->GetXmin();
    const Double_t  xmax  = target.GetXaxis()->GetXmax();
    const Int_t     npar  = (Int_t)fitParam.size();

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

    //
    // setup fit model function
    //

    ReweightPDFFunc modelFunc( coefNames, fitParam, source, sourceCoefs, sourceEval );

    TF1 fitFunc( "EFT", &modelFunc, xmin, xmax, npar );

    for (Int_t i = 0; i < npar; ++i)
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

    std::string fitOption1 = "N Q";
    std::string fitOption2 = "N S E";  // E=hesse error (invert error matrix) and run minos errors

    // use log-likelihood for TH1D but not TProfile
    // Note: do not use "WL" as it does something weird

    bool bLogLike = fitData.IsHist();
    if (bLogLike)
    {
        fitOption1 += " L";
        fitOption2 += " L";
    }

    //
    // Setup minimizer options
    //

    // increase error definition to 2 sigma
    // error definition: the objective function delta that determines parameter error results
    // default is 1 sigma - likelihood: 0.5  chi^2: 1.0
    // both have quadratic dependence, so 2 sigma is 4 times the 1 sigma error definition
    ROOT::Math::MinimizerOptions::SetDefaultErrorDef( bLogLike ? 2 : 4 );

    //
    // first fit with limits
    //

    TFitResultPtr fitStatus;

    fitStatus = fitData.Fit( &fitFunc, fitOption1.c_str() );
    if ((int)fitStatus != 0)
    {
        if (((int)fitStatus < 0) || ((int)fitStatus % 1000 != 0))  // ignore improve (M) errors
            return FitResult( (int)fitStatus );
    }

    //LogMsgInfo( "-*-*-*-*-*-" );

    //
    // second fit with no limits
    //

    // remove parameter limits on non-fixed parameters
    for (Int_t i = 0; i < npar; ++i)
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

    FitResult result;
    {
        result.chi2     = fitStatus->Chi2();
        result.ndf      = fitStatus->Ndf();
        result.prob     = fitStatus->Prob();
        result.chi2_ndf = (result.ndf > 0 ? result.chi2 / result.ndf : 0);

        // get last Fitter object used internally by last Fit call
        TVirtualFitter  *       pFitter      = TVirtualFitter::GetFitter();
        TBackCompFitter *       pBackFitter  = pFitter->InheritsFrom(TBackCompFitter::Class()) ? static_cast<TBackCompFitter *>(pFitter) : nullptr;
        ROOT::Math::Minimizer * pMinimizer   = pBackFitter ? pBackFitter->GetMinimizer() : nullptr;

        for (Int_t i = 0; i < npar; ++i)
        {
            if ((fitIndex >= 0) && (i != fitIndex))
                continue;

            FitResult::Param resultParam;

            resultParam.name  = fitParam[i].name;
            resultParam.value = fitStatus->Parameter(i);
            resultParam.error = fitStatus->ParError(i);

            if (fitStatus->HasMinosError(i))
            {
                resultParam.minosError[0] = fitStatus->LowerError(i);
                resultParam.minosError[1] = fitStatus->UpperError(i);
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

                    std::string sName  = ((fitIndex < 0) ? "min_all_"  : "min_one_") + std::string(obs.name)  + "_" + std::string(fitParam[i].name);
                    std::string sTitle = std::string(obs.title) + ": fit min. for " + std::string(fitParam[i].name) + ((fitIndex < 0) ? " (all)" : " (one)");

                    pGraph->SetName(  sName .c_str() );
                    pGraph->SetTitle( sTitle.c_str() );

                    pGraph->GetXaxis()->SetTitle( fitParam[i].name );
                    pGraph->GetYaxis()->SetTitle( bLogLike ? "Log likelihood" : "#chi^{2}" );
                }
            }

            result.param.push_back( std::move(resultParam) );
        }
    }

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

    typedef std::pair<const char *, int>            NameIndexPair;
    typedef std::map< NameIndexPair, FitResult >    FitResultMap;

    FitResultMap results;

    size_t fitFail = 0;

    for (size_t obsIndex = 0; obsIndex < observables.size(); ++obsIndex)
    {
        const ModelCompare::Observable &    obs         = observables[obsIndex];
        const TH1D *                        pFitTarget  = targetHists[obsIndex];
        const TH1D &                        source      = *sourceHists[obsIndex];
        const ConstTH1DVector               srcCoefs    = ToConstTH1DVector(sourceCoefs[obsIndex]);

        ModelCompare::GoodBadHists goodBad;

        WriteLog( fpLog, "\n******************** %hs ********************", FMT_HS(obs.name) );

        if (pFitTarget->InheritsFrom(TProfile::Class()))
        {
            goodBad = ModelCompare::HistSplitGoodBadBins( pFitTarget,         rawtargetHists[obsIndex] );
            goodBad = ModelCompare::HistSplitGoodBadBins( goodBad.good.get(), rawsourceHists[obsIndex] );
            pFitTarget = goodBad.good.get();
        }

        for (int fitIndex = 0; fitIndex < fitParam.size(); ++fitIndex)
        {
            WriteLog( fpLog, "\n------------------------------------------------------------" );
            WriteLog( fpLog, "Fitting %hs to parameter %hs",
                                FMT_HS(obs.name), FMT_HS(fitParam[fitIndex].name) );
            LogMsgInfo( "------------------------------------------------------------" );

            FitResult fitResult = FitEFTObs( obs, coefNames, *pFitTarget, fitParam,
                                             source, srcCoefs, sourceEval,
                                             fitIndex );

            if (!WriteFitResult( fpLog, fitResult, upOutputFile.get() ))
                ++fitFail;

            // cross-check
            CrossCheckFitResult( coefNames, fitResult, fitParam, *pFitTarget,
                                 source, srcCoefs, sourceEval );

            results[ NameIndexPair(obs.name,fitIndex) ] = std::move(fitResult);
        }

        WriteLog( fpLog, "\n------------------------------------------------------------" );
        WriteLog( fpLog, "Fitting %hs to all parameters",
                            FMT_HS(obs.name) );
        LogMsgInfo( "------------------------------------------------------------") ;

        FitResult fitResult = FitEFTObs( obs, coefNames, *pFitTarget, fitParam,
                                         source, srcCoefs, sourceEval,
                                         -1 );

        if (!WriteFitResult( fpLog, fitResult, upOutputFile.get() ))
            ++fitFail;

        // cross-check
        CrossCheckFitResult( coefNames, fitResult, fitParam, *pFitTarget,
                             source, srcCoefs, sourceEval );

        results[ NameIndexPair(obs.name,-1) ] = std::move(fitResult);

        // TODO: make figures from fit result
    }

    WriteLog( fpLog, "\nSUMMARY" );

    // Write summary
    for (int fitIndex = 0; fitIndex < fitParam.size(); ++fitIndex)
    {
        const FitParam & fpar = fitParam[fitIndex];

        WriteLog( fpLog, "\n---------- %hs ----------", FMT_HS(fpar.name) );

        for (size_t obsIndex = 0; obsIndex < observables.size(); ++obsIndex)
        {
            const ModelCompare::Observable & obs = observables[obsIndex];

            auto itrRes = results.find( NameIndexPair(obs.name,fitIndex) );
            if (itrRes == results.end())
                continue;

            const FitResult & fitOne = itrRes->second;

            itrRes = results.find( NameIndexPair(obs.name,-1) );
            if (itrRes == results.end())
                continue;

            const FitResult & fitAll = itrRes->second;

            auto itrPar = std::find_if( fitOne.param.cbegin(), fitOne.param.cend(), [&](const FitResult::Param & p) { return strcmp(p.name, fpar.name) == 0; } );
            if (itrPar == fitOne.param.cend())
                continue;

            const FitResult::Param & parOne = *itrPar;

            itrPar = std::find_if( fitAll.param.cbegin(), fitAll.param.cend(), [&](const FitResult::Param & p) { return strcmp(p.name, fpar.name) == 0; } );
            if (itrPar == fitAll.param.cend())
                continue;

            const FitResult::Param & parAll = *itrPar;

            auto GetMinos = [](const FitResult::Param & par, double scale = 1E6) -> std::string
            {
                double err = par.error         * scale;
                double lwr = par.minosError[0] * scale;
                double upr = par.minosError[1] * scale;

                if ((lwr == 0.0) && (upr == 0.0)) return std::string();

                std::string sErr = StringFormat( "%.2g", FMT_F(err)  );
                std::string sLwr = StringFormat( "%.2g", FMT_F(-lwr) );
                std::string sUpr = StringFormat( "%.2g", FMT_F(upr)  );

                if ((sLwr == sErr) && (sUpr == sErr)) return std::string();

                return StringFormat( "(%.2g, %.2g)", FMT_F(lwr), FMT_F(upr) );
            };

            std::string minosOne = GetMinos(parOne);
            std::string minosAll = GetMinos(parAll);

            WriteLog( fpLog, "%12hs: %8.2g | %6.2g %16hs | %6.2g || %8.2g | %6.2g %16hs | %6.2g"

                           //"    || %.15E | %.15E | %.15E | %.15E"

                        , FMT_HS(obs.name)
                        , FMT_F(parOne.value*1E6), FMT_F(parOne.error*1E6), FMT_HS(minosOne.c_str()), FMT_F(fitOne.chi2_ndf)
                        , FMT_F(parAll.value*1E6), FMT_F(parAll.error*1E6), FMT_HS(minosAll.c_str()), FMT_F(fitAll.chi2_ndf)

                      //, FMT_F(parOne.value*1E6), FMT_F(parOne.error*1E6)
                      //, FMT_F(parAll.value*1E6), FMT_F(parAll.error*1E6)
                    );
        }
    }

    if (fitFail)
        WriteLog( fpLog, "\n!!! Error: %u fit(s) failed !!!", FMT_U(fitFail) );

    WriteLog( fpLog, "" );

    fclose(fpLog);
}

////////////////////////////////////////////////////////////////////////////////

}  // namespace FitEFT
