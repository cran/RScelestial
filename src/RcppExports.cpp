// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// _scelestial
List _scelestial(DataFrame data, int minK, int maxK);
RcppExport SEXP _RScelestial__scelestial(SEXP dataSEXP, SEXP minKSEXP, SEXP maxKSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< DataFrame >::type data(dataSEXP);
    Rcpp::traits::input_parameter< int >::type minK(minKSEXP);
    Rcpp::traits::input_parameter< int >::type maxK(maxKSEXP);
    rcpp_result_gen = Rcpp::wrap(_scelestial(data, minK, maxK));
    return rcpp_result_gen;
END_RCPP
}
// _synthesis
List _synthesis(int sample, int site, int evolutionSteps, double mutationRate, double advantageIncreaseRatio, double advantageDecreaseRatio, double advantageKeepRatio, double advantageIncreaseStep, double advantageDecreaseStep, double mvRate, double fpRate, double fnRate, int seed);
RcppExport SEXP _RScelestial__synthesis(SEXP sampleSEXP, SEXP siteSEXP, SEXP evolutionStepsSEXP, SEXP mutationRateSEXP, SEXP advantageIncreaseRatioSEXP, SEXP advantageDecreaseRatioSEXP, SEXP advantageKeepRatioSEXP, SEXP advantageIncreaseStepSEXP, SEXP advantageDecreaseStepSEXP, SEXP mvRateSEXP, SEXP fpRateSEXP, SEXP fnRateSEXP, SEXP seedSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type sample(sampleSEXP);
    Rcpp::traits::input_parameter< int >::type site(siteSEXP);
    Rcpp::traits::input_parameter< int >::type evolutionSteps(evolutionStepsSEXP);
    Rcpp::traits::input_parameter< double >::type mutationRate(mutationRateSEXP);
    Rcpp::traits::input_parameter< double >::type advantageIncreaseRatio(advantageIncreaseRatioSEXP);
    Rcpp::traits::input_parameter< double >::type advantageDecreaseRatio(advantageDecreaseRatioSEXP);
    Rcpp::traits::input_parameter< double >::type advantageKeepRatio(advantageKeepRatioSEXP);
    Rcpp::traits::input_parameter< double >::type advantageIncreaseStep(advantageIncreaseStepSEXP);
    Rcpp::traits::input_parameter< double >::type advantageDecreaseStep(advantageDecreaseStepSEXP);
    Rcpp::traits::input_parameter< double >::type mvRate(mvRateSEXP);
    Rcpp::traits::input_parameter< double >::type fpRate(fpRateSEXP);
    Rcpp::traits::input_parameter< double >::type fnRate(fnRateSEXP);
    Rcpp::traits::input_parameter< int >::type seed(seedSEXP);
    rcpp_result_gen = Rcpp::wrap(_synthesis(sample, site, evolutionSteps, mutationRate, advantageIncreaseRatio, advantageDecreaseRatio, advantageKeepRatio, advantageIncreaseStep, advantageDecreaseStep, mvRate, fpRate, fnRate, seed));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_RScelestial__scelestial", (DL_FUNC) &_RScelestial__scelestial, 3},
    {"_RScelestial__synthesis", (DL_FUNC) &_RScelestial__synthesis, 13},
    {NULL, NULL, 0}
};

RcppExport void R_init_RScelestial(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
