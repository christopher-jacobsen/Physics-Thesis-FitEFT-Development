//
//  FitEFT.h
//  FitEFT
//
//  Created by Christopher Jacobsen on 23/09/15.
//  Copyright (c) 2015 Christopher Jacobsen. All rights reserved.
//

#ifndef FIT_EFT_H
#define FIT_EFT_H

#include "RootUtil.h"
#include "ModelCompare.h"
#include "ReweightEFT.h"

////////////////////////////////////////////////////////////////////////////////

namespace FitEFT
{

////////////////////////////////////////////////////////////////////////////////

void LoadTupleData( const ModelCompare::ObservableVector & observables, const ModelCompare::ModelFile & model,
                    RootUtil:: TupleVector & modelData,
                    const char * cacheFileName );

void LoadTupleData( const ModelCompare::ObservableVector & observables, const ModelCompare::ModelFileVector & models,
                    std::vector<RootUtil::TupleVector> & allData,
                    const char * cacheFileName );

void MergeZeroBins( TH1D & hist, bool bOnlyInner = true );

////////////////////////////////////////////////////////////////////////////////

struct FitParam
{
    const char * name;
    double       initValue;
    double       minValue;
    double       maxValue;
};

typedef std::vector<FitParam> FitParamVector;

////////////////////////////////////////////////////////////////////////////////

enum class FitKind
{
    Undefined,
    Binned,
    Binned_ShapeOnly,
    Unbinned
};

////////////////////////////////////////////////////////////////////////////////

void FitEFT( const char * outputFileName,
             const ModelCompare::ObservableVector & observables,
             const RootUtil::CStringVector & coefNames,
             const ModelCompare::ModelFile & targetFile, const FitParamVector & fitParam,
             const ModelCompare::ModelFile & sourceFile, const ReweightEFT::ParamVector & sourceParam,
             double luminosity,
             const std::vector<FitKind> & fitKinds, bool bFitAll, bool bCreateScanGraph,
             const char * cacheFileName );

////////////////////////////////////////////////////////////////////////////////

}  // namespace FitEFT

#endif // FIT_EFT_H
