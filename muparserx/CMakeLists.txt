# ==============================================================================
#
#                                 r  i  t  a
#
#            An environment for Modelling and Numerical Simulation
#
# ==============================================================================
#
#   Copyright (C) 2021 - 2022  Rachid Touzani
#
#   This file is part of rita.
#
#   rita is free software: you can redistribute it and/or modify
#   it under the terms of the GNU General Public License as published by
#   the Free Software Foundation, either version 3 of the License, or
#   (at your option) any later version.
#
#   rita is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.
#
# ==============================================================================

target_sources (${PROJECT_NAME} PRIVATE
                mpError.cpp
                mpOprtPostfixCommon.cpp
                mpFuncCmplx.cpp
                mpPackageCmplx.cpp
                mpFuncCommon.cpp
                mpPackageCommon.cpp
                mpFuncMatrix.cpp
                mpPackageMatrix.cpp
                mpFuncNonCmplx.cpp
                mpPackageNonCmplx.cpp
                mpFuncStr.cpp
                mpPackageStr.cpp
                mpICallback.cpp
                mpPackageUnit.cpp
                mpIOprt.cpp
                mpParser.cpp
                mpIPackage.cpp
                mpParserBase.cpp
                mpIToken.cpp
                mpParserMessageProvider.cpp
                mpIValReader.cpp
                mpRPN.cpp
                mpIValue.cpp
                mpScriptTokens.cpp
                mpIfThenElse.cpp
                mpTest.cpp
                mpOprtBinAssign.cpp
                mpTokenReader.cpp
                mpOprtBinCommon.cpp
                mpValReader.cpp
                mpOprtCmplx.cpp
                mpValue.cpp
                mpOprtIndex.cpp
                mpValueCache.cpp
                mpOprtMatrix.cpp
                mpVariable.cpp
                mpOprtNonCmplx.cpp
               )
