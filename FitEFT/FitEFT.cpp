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
#include <TF1.h>
#include <TFitResult.h>
#include <Math/MinimizerOptions.h>

using namespace RootUtil;

////////////////////////////////////////////////////////////////////////////////

namespace FitEFT
{

struct FitResult
{
    struct Param
    {
        const char *    name;
        double          value;
        double          error;
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
FitResult FitEFTObs( const ModelCompare::Observable & obs, const CStringVector & coefNames,
                     const TH1D & target, const FitParamVector & fitParam,
                     const TH1D & source, const ConstTH1DVector & sourceCoefs, const std::vector<double> & sourceEval,
                     int fitIndex = -1 )   // -1 = fit all, otherwise index of fit parameter to fit, keeping all others fixed
{
// "L" : error independent of scale
//   ocWWW =  1.407493792E-06 ± 3.929595539E-07
//     ocW = -3.082096616E-07 ± 3.933420127E-07
//     ocB = -7.768707211E-05 ± 3.818847769E-05
// all:
//   ocWWW = -1.133606645E-07 ± 9.302806678E-07
//     ocW = -1.225535299E-06 ± 4.127463750E-07
//     ocB = -8.837049483E-05 ± 2.475648113E-05

// "WL" : error depends on scale for some parameters
//   ocWWW =  1.407496157E-06 ± 1.535264034E-xx
//     ocW = -3.082091173E-07 ± 3.105637463E-07
//     ocB = -7.768706762E-05 ± 1.470158572E-05
// all:
//   ocWWW = -1.133630765E-07 ± 5.948353820E+xx
//     ocW = -1.225537852E-06 ± 1.840232309E+xx
//     ocB = -8.837058852E-05 ± 6.133336659E+xx

  // PTZ:
  //const Double_t  scale = 1;      // errors =    E-13,    E-7,    E-5 /    E-14    E-15    E-13 (error)
  //const Double_t  scale = 1E-3;   // errors = 1.5E-10, 3.1E-7, 1.5E-5 / 4.9E-4  1.5E-4  5.1E-2  (warn)
  //const Double_t  scale = 1E-6;   // errors = 1.5E-7,  3.1E-7, 1.5E-5 / 1.6E-5  5.1E-6  1.7E-3  (warn)
  //const Double_t  scale = 1E-9;   // errors = 1.5E-4,  3.1E-7, 1.5E-5 / 5.9E-3  1.8E-3  6.1E-1
  //const Double_t  scale = 1E-10;  // errors =    E-3,     E-7,    E-5 /    E-2     E-2     E+0
  //const Double_t  scale = 1E-11;  // errors =    E-2,     E-7,    E-5 /    E-1     E-1     E+1
  //const Double_t  scale = 1E-12;  // errors = 1.5E-1,  3.1E-7, 1.5E-5 / 5.9E+0  1.8E+0  6.1E+2
  //const Double_t  scale = 1E-16;  // errors =    E-3,     E-7,    E-5 /    E+4     E+4     E+6

  // O1(cWWW)
//   ocWWW =  6.273085662E-09 ± 2.436393358E-17
//     ocW =  6.982349016E-07 ± 4.581396056E-17
//     ocB = -9.385725459E-06 ± 3.284975815E-15
//
//   ocWWW = -3.252253448E-10 ± 2.755932386E-17
//     ocW =  6.960137960E-07 ± 4.635048099E-17
//     ocB = -1.194507236E-06 ± 3.736071058E-15
    const Double_t  scale = 1;      // errors = 2.4E-17,  4.6E-17, 3.3E-17 / 4.9E-18  4.6E-17  3.7E-15
  //const Double_t  scale = 1E-6;   // errors = 2.4E-17,  4.6E-17, (Fail)  / 2.8E-17  4.6E-17  3.8E-15

    const Double_t  xmin  = target.GetXaxis()->GetXmin();
    const Double_t  xmax  = target.GetXaxis()->GetXmax();
    const Int_t     npar  = (Int_t)fitParam.size();

    /////

    // define and setup fit lambda function

    size_t                      rejectCount(0);
    std::vector<Double_t>       lastPar(npar);
    ConstTH1DUniquePtr          upLastSourceReweight;
    ReweightEFT::ParamVector    targetParam(npar);

    for (Int_t i = 0; i < npar; ++i)
    {
        targetParam[i].name  = fitParam[i].name;
        targetParam[i].value = fitParam[i].initValue;
    }

    auto EFTFunc = [&](const Double_t * x, const Double_t * par) -> Double_t
    {
        if (!upLastSourceReweight ||
            (memcmp( lastPar.data(), par, npar * sizeof(*par)) != 0))
        {
            for (Int_t i = 0; i < npar; ++i)
                targetParam[i].value = par[i] * scale;

            std::vector<double> targetEval;
            CalcEvalVector( coefNames, targetParam, targetEval );

            upLastSourceReweight.reset(
                ReweightEFT::ReweightHist( source, sourceCoefs, sourceEval, targetEval,
                                           "RWFitHist", "RWFitHist" ) );

            memcpy( lastPar.data(), par, npar * sizeof(*par) );
        }

        const TH1D * pHist = upLastSourceReweight.get();

        Double_t xVal    = x[0];
        Int_t    bin     = pHist->FindFixBin( xVal );
        Double_t nEff    = GetHistBinEffectiveEntries( *pHist, bin );
        Double_t content = pHist->GetBinContent( bin );

        //LogMsgInfo( "x=%f b=%i n=%f c=%f", FMT_F(xVal), FMT_I(bin), FMT_F(nEff), FMT_F(content) );

        // reduce NDF
        if (nEff <= 0)
        {
            // LogMsgInfo( ">>>> Rejecting bin %i <<<<", bin );
            ++rejectCount;
            TF1::RejectPoint();
            return 0;
        }

        return content;
    };

    //////

    // setup fit function

    TF1 fitFunc( "EFT", EFTFunc, xmin, xmax, npar );

    for (Int_t i = 0; i < npar; ++i)
    {
        fitFunc.SetParName(   i, fitParam[i].name );
        fitFunc.SetParameter( i, fitParam[i].initValue / scale );
        fitFunc.SetParLimits( i, fitParam[i].minValue / scale, fitParam[i].maxValue / scale );
	
        // set initial step (helps fit converge)
        fitFunc.SetParError(  i, (fitParam[i].maxValue - fitParam[i].minValue) / 1E6 / scale );

        if ((fitIndex >= 0) && (i != fitIndex))
            fitFunc.FixParameter( i, fitParam[i].initValue / scale );
    }

    TH1DUniquePtr pFitHist{ (TH1D *)target.Clone() };     // clone hist as Fit is not const

    // setup fit options
    // do not use M - it produces too many warnings of:
    // "FUNCTION VALUE DOES NOT SEEM TO DEPEND ON ANY OF THE 1 VARIABLE PARAMETERS. VERIFY THAT STEP SIZES ARE BIG ENOUGH AND CHECK FCN LOGIC."

    std::string fitOption1 = "N Q";
    std::string fitOption2 = "N S";

    // use log-likelihood for TH1D but not TProfile
    // Note: do not use "WL" as it does something weird

    if (!target.InheritsFrom(TProfile::Class()))
    {
        fitOption1 += " L";
        fitOption2 += " L";
    }

    TFitResultPtr fitStatus;

    //
    // first fit with limits
    //

    fitStatus = pFitHist->Fit( &fitFunc, fitOption1.c_str() );
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

    rejectCount = 0;

    // re-run fit using result from previous fit as initial values
    // including parameter errors that are used to setup fit steps
    fitStatus = pFitHist->Fit( &fitFunc, fitOption2.c_str() );
    if ((int)fitStatus != 0)
    {
        if (((int)fitStatus < 0) || ((int)fitStatus % 1000 != 0))   // ignore improve (M) errors
            return FitResult( (int)fitStatus );
    }

    fitStatus->Print();

    LogMsgInfo( "\nReject Count: %u", FMT_U(rejectCount) );

    // fill in result

    FitResult result;
    {
        result.chi2     = fitFunc.GetChisquare();
        result.ndf      = fitFunc.GetNDF();
        result.prob     = fitFunc.GetProb();
        result.chi2_ndf = (result.ndf > 0 ? result.chi2 / result.ndf : 0);

        for (Int_t i = 0; i < npar; ++i)
        {
            if ((fitIndex >= 0) && (i != fitIndex))
                continue;

            FitResult::Param resultParam;

            resultParam.name  = fitParam[i].name;
            resultParam.value = fitFunc.GetParameter(i) * scale;
            resultParam.error = fitFunc.GetParError(i)  * scale;

            result.param.push_back(resultParam);
        }
    }

    // cross-check chi2
    {
        Double_t chi2 = target.Chisquare( &fitFunc );

        LogMsgInfo( "Cross-check chi2: %g", FMT_F(chi2) );
    }

    LogMsgHistBinCounts( target );

    LogMsgInfo( "" );  // empty line

    return result;
}

////////////////////////////////////////////////////////////////////////////////
TH1D * CreateReweightFitResult( const CStringVector & coefNames,
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

        for (const auto & t : targetParam)
        {
            LogMsgInfo( "> %hs = %g", FMT_HS(t.name), FMT_F(t.value) );
        }
    }

    // construct target evaluation
    std::vector<double> targetEval;
    CalcEvalVector( coefNames, targetParam, targetEval );

    TH1D * pHist = ReweightEFT::ReweightHist( source, sourceCoefs, sourceEval, targetEval, "RWFitHist", "RWFitHist" );

    return pHist;
}

////////////////////////////////////////////////////////////////////////////////
void CrossCheckFitResult( const CStringVector & coefNames,
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

    upFitHist.reset( ConvertTProfileToTH1D( upFitHist.get(), false ) );

    // zero errors on reweightHist (make perfect)
    {
        const Int_t nSize = upFitHist->GetSize();
        for (Int_t bin = 0; bin < nSize; ++bin)
        {
            upFitHist->SetBinError(bin, 0.0);
        }
        upFitHist->ResetStats();
        //LogMsgHistDump(*upFitHist);
    }

    TH1DUniquePtr upChiHist{ (TH1D *)target.Clone("ChiHist") };

    // zero bins that have zero error in chi hist
    {
        const Int_t nSize = upChiHist->GetSize();
        for (Int_t bin = 0; bin < nSize; ++bin)
        {
            Double_t error = upChiHist->GetBinError(bin);
            if (error == 0)
            {
                ZeroHistBin( *upChiHist, bin );
                ZeroHistBin( *upFitHist, bin );
            }
        }
        upChiHist->ResetStats();
    }

    LogMsgHistEffectiveEntries( *upChiHist );
    LogMsgHistEffectiveEntries( *upFitHist );

    LogMsgHistBinCounts( *upChiHist );
    LogMsgHistBinCounts( *upFitHist );

    LogMsgInfo( "" );

    RootUtil::Chi2Result chi2;
    chi2.Chi2Test( *upChiHist, *upFitHist );

    LogMsgInfo( chi2.Label().c_str() );
}

////////////////////////////////////////////////////////////////////////////////
void WriteLog( FILE * fpLog, const char * format, ... )
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
bool WriteFitResult( FILE * fpLog, const FitResult & result )
{
    if (result.status != 0)
    {
        WriteLog( fpLog, "!!! Fit Failed: status = %i", FMT_I(result.status) );
        return false;
    }

    WriteLog( fpLog, "------------------------------------------------------------" );

    WriteLog( fpLog, "chi2 / ndf = %g / %g = %g", FMT_F(result.chi2), FMT_F(result.ndf), FMT_F(result.chi2_ndf) );
    WriteLog( fpLog, "prob = %g\n", FMT_F(result.prob) );

    for (const FitResult::Param & p : result.param)
    {
        WriteLog( fpLog, "%8hs = %.9E ± %.9E", FMT_HS(p.name), FMT_F(p.value), FMT_F(p.error) );
    }

    WriteLog( fpLog, "------------------------------------------------------------" );

    return true;
}

////////////////////////////////////////////////////////////////////////////////
void FitEFT( const char * outputFileName,
             const ModelCompare::ObservableVector & observables,
             const CStringVector & coefNames,
             const ModelCompare::ModelFile & targetFile, const FitParamVector & fitParam,
             const ModelCompare::ModelFile & sourceFile, const ReweightEFT::ParamVector & sourceParam )
{
    // disable automatic histogram addition to current directory
    TH1::AddDirectory(kFALSE);
    // enable automatic sumw2 for every histogram
    TH1::SetDefaultSumw2(kTRUE);

    //ROOT::Math::MinimizerOptions::SetDefaultMinimizer("Minuit2");

    //ROOT::Math::MinimizerOptions minOpt;
    //minOpt.Print();

    //ROOT::Math::MinimizerOptions::PrintDefault("Minuit2");

    //return;

    // ------

    std::string rootFileName = std::string(outputFileName) + ".root";
    std::string logFileName  = std::string(outputFileName) + ".txt";

    LogMsgInfo( "Root output file: %hs", FMT_HS(rootFileName.c_str()) );
    std::unique_ptr<TFile> upOutputFile( new TFile( rootFileName.c_str(), "RECREATE" ) );
    if (upOutputFile->IsZombie() || !upOutputFile->IsOpen())    // IsZombie is true if constructor failed
    {
        LogMsgError( "Failed to create output file (%hs).", FMT_HS(rootFileName.c_str()) );
        ThrowError( std::invalid_argument( rootFileName ) );
    }

    LogMsgInfo( "Log output file: %hs", FMT_HS(logFileName.c_str()) );
    FILE * fpLog = fopen( logFileName.c_str(), "wt" );

    // load input histograms

    TH1DVector              targetData;     // targetData[observable]
    TH1DVector              sourceData;     // sourceData[observable]
    std::vector<TH1DVector> sourceCoefs;    // sourceCoefs[observable][coefficient]
    std::vector<double>     sourceEval;     // sourceEval[coefficient]
    TH1DVector              rawTargetData;  // rawTargetData[observable]
    TH1DVector              rawSourceData;  // rawSourceData[observable]

    ReweightEFT::LoadReweightFiles( observables, coefNames, targetFile, sourceFile, sourceParam,    // inputs
                                    targetData, sourceData, sourceCoefs, sourceEval,                // outputs
                                    rawTargetData, rawSourceData );

    // write input hists

    WriteHists( upOutputFile.get(), targetData );   // output file takes ownership of histograms
    WriteHists( upOutputFile.get(), sourceData );   // output file takes ownership of histograms

    // uncomment to save coefficent hists to file
    //  for (const auto & coefData : sourceCoefs)
    //      WriteHists( upOutputFile.get(), coefData );  // output file takes ownership of histograms

    LogMsgInfo("");
    LogMsgHistBinCounts( ToConstTH1DVector(targetData) );
    LogMsgInfo("");
    LogMsgHistBinCounts( ToConstTH1DVector(sourceData) );
    LogMsgInfo("");
    LogMsgHistBinCounts( ToConstTH1DVector(targetData), ToConstTH1DVector(sourceData) );
    LogMsgInfo("");

    /*
    // replace targetData errors with sqrt(N) errors - Asimov samples
    // The current errors on targetData, are errors on the expectation value.
    // To convert this to pseudo-data, we keep the expectation value, but replace
    // the errors with sqrt(binContent).

    LogMsgInfo("------------ Before Asimov ------------");
    targetData[0]->Print("all");
    LogMsgHistStats(*targetData[0]);

    for (TH1D * pHist : targetData)
    {
        pHist->Sumw2(kFALSE);
        pHist->Sumw2(kTRUE);
        pHist->ResetStats();
    }

    LogMsgInfo("------------ After Asimov ------------");
    LogMsgHistStats(*targetData[0]);
    targetData[0]->Print("all");
    */

    // fit each source hist to the target, varying only one parameter

    size_t fitFail = 0;

    for (size_t obsIndex = 0; obsIndex < observables.size(); ++obsIndex)
    {
        const ModelCompare::Observable &    obs         = observables[obsIndex];
        const TH1D *                        pFitTarget  = targetData[obsIndex];
        const TH1D &                        source      = *sourceData[obsIndex];
        const ConstTH1DVector               srcCoefs    = ToConstTH1DVector(sourceCoefs[obsIndex]);

        ModelCompare::GoodBadHists goodBad;

        WriteLog( fpLog, "\n******************** %hs ********************", FMT_HS(obs.name) );

        if (pFitTarget->InheritsFrom(TProfile::Class()))
        {
            goodBad = ModelCompare::HistSplitGoodBadBins( pFitTarget, rawTargetData[obsIndex] );
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

            if (!WriteFitResult( fpLog, fitResult ))
                ++fitFail;

            // cross-check
            CrossCheckFitResult( coefNames, fitResult, fitParam, *pFitTarget,
                                 source, srcCoefs, sourceEval );
        }

        WriteLog( fpLog, "\n------------------------------------------------------------" );
        WriteLog( fpLog, "Fitting %hs to all parameters",
                            FMT_HS(obs.name) );
        LogMsgInfo( "------------------------------------------------------------") ;

        FitResult fitResult = FitEFTObs( obs, coefNames, *pFitTarget, fitParam,
                                         source, srcCoefs, sourceEval,
                                         -1 );

        if (!WriteFitResult( fpLog, fitResult ))
            ++fitFail;

        // cross-check
        CrossCheckFitResult( coefNames, fitResult, fitParam, *pFitTarget,
                             source, srcCoefs, sourceEval );

        // TODO: make figures from fit result
    }

    if (fitFail)
        WriteLog( fpLog, "\n!!! Error: %u fit(s) failed !!!", FMT_U(fitFail) );

    fclose(fpLog);
}

////////////////////////////////////////////////////////////////////////////////

}  // namespace FitEFT
