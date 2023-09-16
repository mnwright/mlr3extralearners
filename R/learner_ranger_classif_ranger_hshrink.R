#' @title Classification Rf With Hierarchical Shrinkage Learner
#' @author mnwright
#' @name mlr_learners_classif.ranger_hshrink
#'
#' @description
#' FIXME: BRIEF DESCRIPTION OF THE LEARNER.
#' Calls [ranger::ranger()] from FIXME: (CRAN VS NO CRAN): \CRANpkg{ranger} | 'ranger'.
#'
#' @section Initial parameter values:
#' FIXME: DEVIATIONS FROM UPSTREAM PARAMETERS. DELETE IF NOT APPLICABLE.
#'
#' @section Custom mlr3 parameters:
#' FIXME: DEVIATIONS FROM UPSTREAM DEFAULTS. DELETE IF NOT APPLICABLE.
#'
#' @templateVar id classif.ranger_hshrink
#' @template learner
#'
#' @references
#' `r format_bib(FIXME: ONE OR MORE REFERENCES FROM bibentries.R)`
#'
#' @template seealso_learner
#' @template example
#' @export
LearnerClassifRangerHshrink = R6Class("LearnerClassifRangerHshrink",
  inherit = LearnerClassif,
  public = list(
    #' @description
    #' Creates a new instance of this [R6][R6::R6Class] class.
    initialize = function() {
      param_set = ps(
        # FIXME: only works if predict_type == "se". How to set dependency?
        alpha                        = p_dbl(default = 0.5, tags = "train"),
        always.split.variables       = p_uty(tags = "train"),
        class.weights                = p_uty(default = NULL, tags = "train"),
        holdout                      = p_lgl(default = FALSE, tags = "train"),
        importance                   = p_fct(c("none", "impurity", "impurity_corrected", "permutation"), tags = "train"),
        keep.inbag                   = p_lgl(default = FALSE, tags = "train"),
        max.depth                    = p_int(default = NULL, lower = 0L, special_vals = list(NULL), tags = "train"),
        min.bucket                   = p_int(1L, default = 1L, tags = "train"),
        min.node.size                = p_int(1L, default = NULL, special_vals = list(NULL), tags = "train"),
        min.prop                     = p_dbl(default = 0.1, tags = "train"),
        minprop                      = p_dbl(default = 0.1, tags = "train"),
        mtry                         = p_int(lower = 1L, special_vals = list(NULL), tags = "train"),
        mtry.ratio                   = p_dbl(lower = 0, upper = 1, tags = "train"),
        num.random.splits            = p_int(1L, default = 1L, tags = "train"),
        num.threads                  = p_int(1L, default = 1L, tags = c("train", "predict", "threads")),
        num.trees                    = p_int(1L, default = 500L, tags = c("train", "predict", "hotstart")),
        oob.error                    = p_lgl(default = TRUE, tags = "train"),
        regularization.factor        = p_uty(default = 1, tags = "train"),
        regularization.usedepth      = p_lgl(default = FALSE, tags = "train"),
        replace                      = p_lgl(default = TRUE, tags = "train"),
        respect.unordered.factors    = p_fct(c("ignore", "order", "partition"), default = "ignore", tags = "train"),
        sample.fraction              = p_dbl(0L, 1L, tags = "train"),
        save.memory                  = p_lgl(default = FALSE, tags = "train"),
        scale.permutation.importance = p_lgl(default = FALSE, tags = "train"),
        se.method                    = p_fct(c("jack", "infjack"), default = "infjack", tags = "predict"),
        seed                         = p_int(default = NULL, special_vals = list(NULL), tags = c("train", "predict")),
        split.select.weights         = p_uty(default = NULL, tags = "train"),
        splitrule                    = p_fct(c("gini", "extratrees", "hellinger"), default = "gini", tags = "train"),
        verbose                      = p_lgl(default = TRUE, tags = c("train", "predict")),
        write.forest                 = p_lgl(default = TRUE, tags = "train"),
        lambda                       = p_uty(default = 0, tags = "train"),
        measure                      = p_uty(default = NULL, tags = "train")
      )

      param_set$values = list(num.threads = 1L)

      super$initialize(
        id = "classif.ranger_hshrink",
        packages = "ranger",
        feature_types = c("logical", "integer", "numeric", "character", "factor", "ordered"),
        predict_types = c("prob"),
        param_set = param_set,
        properties = c("multiclass", "twoclass", "oob_error"),
        man = "mlr3extralearners::mlr_learners_classif.ranger_hshrink",
        label = "ranger with hshrink"
      )
    }
  ),
  private = list(
    .train = function(task) {
      pv = self$param_set$get_values(tags = "train")
      pv = convert_ratio(pv, "mtry", "mtry.ratio", length(task$feature_names))

      if (!is.null(pv$lambda)) {
        lambda = pv$lambda
        pv$lambda = NULL
        pv$node.stats = TRUE
        pv$keep.inbag = TRUE
      } else {
        lambda = 0
      }

      if (!is.null(pv$measure)) {
        if (inherits(pv$measure, "Measure")) {
          measure = pv$measure
        } else if (is.character(pv$measure)) {
          measure = msr(pv$measure)
        } else {
          stop("Unknown measure type.")
        }
        pv$measure = NULL
      }

      rf = invoke(ranger::ranger,
             dependent.variable.name = task$target_names,
             data = task$data(),
             probability = self$predict_type == "prob",
             case.weights = task$weights$weight,
             .args = pv
      )

      if (is.numeric(lambda) & length(lambda) == 1) {
        if (lambda > 0) {
          ranger::hshrink(rf, lambda = lambda)
        }
        rf
      } else if (is.numeric(lambda) & is.vector(lambda)) {
        if (is.null(measure)) {
          stop("Need measure for optimizing lambda.")
        }

        res = sapply(lambda, function(ll) {
          rf_shrinked = rlang::duplicate(rf)
          ranger::hshrink(rf_shrinked, lambda = ll)

          # Re-calculate OOB error
          oob_idx = ifelse(simplify2array(rf$inbag.counts) == 0, TRUE, NA)
          preds = predict(rf_shrinked, task$data(), predict.all = TRUE)$predictions

          for (i in 1:dim(preds)[2]) {
            preds[, i, ] = oob_idx * preds[, i, ]
          }
          oob_pred = apply(preds, 1:2, mean, na.rm = TRUE)
          oob_pred[oob_pred < 0] = 0

          measure$fun(task$truth(), oob_pred)
        })

        if (measure$minimize) {
          best_lambda = lambda[which.min(res)]
        } else {
          best_lambda = lambda[which.max(res)]
        }

        #message("best lambda: ", best_lambda)

        # Return best shrinked model
        rf_shrinked = rlang::duplicate(rf)
        ranger::hshrink(rf_shrinked, lambda = best_lambda)
        rf_shrinked
      } else {
        stop("Unknown value for lambda.")
      }

    },
    .predict = function(task) {
      pv = self$param_set$get_values(tags = "predict")
      newdata = ordered_features(task, self)

      prediction = invoke(predict,
                          self$model,
                          data = newdata,
                          predict.type = "response", .args = pv
      )

      if (self$predict_type == "response") {
        list(response = prediction$predictions)
      } else {
        list(prob = prediction$predictions)
      }
    }
  )
)

.extralrns_dict$add("classif.ranger_hshrink", LearnerClassifRangerHshrink)
