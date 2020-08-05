/***************************************************************************
* Copyright (c) Wolf Vollprecht, Johan Mabille and Sylvain Corlay          *
* Copyright (c) QuantStack                                                 *
*                                                                          *
* Distributed under the terms of the BSD 3-Clause License.                 *
*                                                                          *
* The full license is in the file LICENSE, distributed with this software. *
****************************************************************************/

#ifndef XBLAS_CONFIG_CLING_HPP
#define XBLAS_CONFIG_CLING_HPP

#pragma cling add_library_path("/usr/local/miniconda3/envs/metapop/lib64")
#pragma cling add_library_path("/usr/local/miniconda3/envs/metapop/lib32")
#pragma cling add_library_path("/usr/local/miniconda3/envs/metapop/lib")

#ifndef XTENSOR_USE_FLENS_BLAS

#define HAVE_CBLAS 1

#pragma cling load("libblas")
#pragma cling load("liblapack")

#endif

#endif
