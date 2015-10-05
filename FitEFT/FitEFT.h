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

struct FitParam
{
    const char * name;
    double       initValue;
    double       minValue;
    double       maxValue;
};

typedef std::vector<FitParam> FitParamVector;

////////////////////////////////////////////////////////////////////////////////

void FitEFT( const char * outputFileName,
             const ModelCompare::ObservableVector & observables,
             const RootUtil::CStringVector & coefNames,
             const ModelCompare::ModelFile & targetFile, const FitParamVector & fitParam,
             const ModelCompare::ModelFile & sourceFile, const ReweightEFT::ParamVector & sourceParam,
             double luminosity,
             const char * cacheFileName );

////////////////////////////////////////////////////////////////////////////////

}  // namespace FitEFT

#endif // FIT_EFT_H
