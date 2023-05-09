#' Applies a baseflow estimation filter using Stan
#'
#' Missing values are allowed in the response but not in the covariates.
#'
#' @param formula A formula as in lm()
#' @param data A long data frame containing the dates, covariates and the response variable.
#' @param method A string specifying the type of filter to use.
#' @param passes The number of passes of the filter. Usually 3.
#' @param time_method Defines the type of time series model to capture temporal dependence.
#' @param iter Number of iterations
#' @param warmup Warm up samples
#' @param chains Number of chains
#' @param refresh Sampler refreshing rate
#' @param seed (optional) A seed for reproducibility
#' @param bfi A baseflow index used in the Eckhardt's approach
#' @return A list with the model fit
#' @details Missing values are not allowed in the covariates and they must be imputed before using the function.
#' @return It returns a stanfit object. It includes the formula used to fit the model.
#' @export
#' @importFrom dplyr mutate %>% distinct left_join case_when
#' @importFrom plyr .
#' @importFrom rstan stan
#' @importFrom stats dist df
#' @author Edgar Santos-Fernandez
#' @examples
#'\dontrun{
#'#options(mc.cores = parallel::detectCores())
#'# Lyne and Hollick filter with 3 passes
#' fit_lh <- baseflow(formula = y ~ prec,
#'                    data = df,
#'                    method = 'lh',
#'                   passes = 3)
#'
#' #Eckhardt approach
#' fit_eck <- baseflow(formula = y ~ prec,
#'                     data = df,
#'                     method = 'eck')
#'}
#'

baseflow <- function(formula = formula,
                     data = data,
                     method = 'lh',
                     passes = 3,
                     bfi = 0.8,
                     time_method = time_method, # list("ar", "date")
                     iter = 3000,
                     warmup = 1500,
                     chains = 3,
                     refresh = max(iter/100, 1),
                     seed = seed
){





  time_method = list("ar", "date") # NB


  data_com <-  'data {
      int<lower=1> T;
      int<lower = 0> N_y_obs; // number observed values
      int<lower = 0> N_y_mis; // number missing values
      int<lower = 1> i_y_obs[N_y_obs] ;  //[N_y_obs,T]
      int<lower = 1> i_y_mis[N_y_mis] ;  // N_y_mis,T]
      vector<lower=0> [N_y_obs] y_obs;
      int<lower=1> passes;
      vector<lower=0> [passes] ii;
      int<lower=1> K;
      matrix[T,K] X ;
      real <lower=0, upper = 1> bfi;
  }'


  param_com <- '
  parameters {
      vector[N_y_mis] y_mis;//declaring the missing y
      real <lower=0, upper = 1> alpha; // NB
      real<lower=0> sigma;
      vector[K] beta;
  //    real <lower=-1, upper = 1> phi; // NB
  '


  param_phi_ar <- '
    real <lower=-1, upper = 1> phi; // NB
  '

  tparam_lh <- '
  transformed parameters {
    vector[T] y; //vector<lower=0> [T] y;
    matrix[T,passes] f; // filter
    matrix[T,passes] b; // vector<lower=0> [T] b; base
    matrix[T,passes+1] y0; // updated y
    vector<lower=0>[T] mu;
    vector[T] epsilon; // error term
      y[i_y_obs] = y_obs;
      y[i_y_mis] = y_mis;

      y0[,1] = y;

  mu[1] = X[1] * beta;
  epsilon[1] = y[1] - mu[1];

  for (t in 2:T){
      mu[t] = X[t] * beta;
      epsilon[t] = y[t] - mu[t];
      mu[t] = mu[t] + phi * epsilon[t-1]; //
  }
      for (i in 1:passes){
          if(i == 1){ y0[1,i] = y[1];}

          // odd passes (forward)
          if(ii[i] == 1){;
            f[1,i] = 0; //y0[1,i] NB
            b[1,i] = y0[1,i];

            for (t in 2:T){
              f[t,i] = alpha * f[t-1,i] + 0.5 * (1 + alpha) * (y0[t,i] - y0[t-1,i]);
              if(f[t,i] > 0){
                    b[t,i] = y0[t,i] - f[t,i]   ;
                  }
                else{
                      b[t,i] = y0[t,i] ;
                  }
            }
           y0[,i+1] = b[,i];

          }



          // even passes (backward)
          if(ii[i] == 2){;
            f[T,i] = 0; //y0[T,i]
            b[T,i] = y0[T,i];
            for (t in 1:T-1){
              f[T-t,i] = alpha * f[T-t+1,i] + 0.5 * (1 + alpha) * (y0[T-t,i] - y0[T-t+1,i]);
              if(f[T-t,i] > 0){
                    b[T-t,i] = y0[T-t,i] - f[T-t,i]   ;
                  }
                else{
                      b[T-t,i] = y0[T-t,i] ;
                  }
            }
            y0[,i+1] = b[,i];
          }

      }
  }

   '

  tparam_eck <- '
  transformed parameters {
    vector[T] y; //vector<lower=0> [T] y;
    matrix[T,passes] f; // filter
    matrix[T,passes] b; // vector<lower=0> [T] b; base
    matrix[T,passes+1] y0; // updated y
    vector<lower=0>[T] mu;
    vector[T] epsilon; // error term
      y[i_y_obs] = y_obs;
      y[i_y_mis] = y_mis;

      y0[,1] = y;

  mu[1] = X[1] * beta;
  epsilon[1] = y[1] - mu[1];

  for (t in 2:T){
      mu[t] = X[t] * beta;
      epsilon[t] = y[t] - mu[t];
      mu[t] = mu[t] + phi * epsilon[t-1]; //
  }
      for (i in 1:passes){
          if(i == 1){ y0[1,i] = y[1];}

          // odd passes (forward)
          if(ii[i] == 1){;
            f[1,i] = 0; //y0[1,i] NB
            b[1,i] = y0[1,i];

            for (t in 2:T){
              f[t,i] = ((1 - bfi) * alpha * f[t-1,i] + (1 - alpha) * bfi * (y0[t,i]) )/ (1 - alpha * bfi) ; //
              if(f[t,i] <=  y0[t,i]){
                    b[t,i] = f[t,i]   ;
                  }
                else{
                      b[t,i] = y0[t,i] ;
                  }
            }
           y0[,i+1] = b[,i];

          }



          // even passes (backward)
          if(ii[i] == 2){;
            f[T,i] = 0; //y0[T,i]
            b[T,i] = y0[T,i];
            for (t in 1:T-1){
              //f[T-t,i] = alpha * f[T-t+1,i] + 0.5 * (1 + alpha) * (y0[T-t,i] - y0[T-t+1,i]);
              f[T-t,i] = ((1 - bfi) * alpha * f[T-t+1,i] + (1 - alpha) * bfi * (y0[T-t+1,i]) ) / (1 - alpha * bfi) ;
              if(f[T-t,i] <=  y0[T-t,i]){
                    b[T-t,i] = f[T-t,i]   ;
                  }
                else{
                      b[T-t,i] = y0[T-t,i] ;
                  }
            }
            y0[,i+1] = b[,i];
          }

      }
  }
  '


  model_com <- '
   model {
      target += normal_lpdf( y | mu, sigma); //y ~ normal(mu, sigma);
      alpha ~ beta(92.5, 7.5); // beta(9, 1) beta(90, 10)
      //mu ~ normal(20,10);
      sigma ~ uniform(0.001,50);
      phi ~ normal(0.5,0.3); // uniform(-1, 1) or can use phi ~ normal(0.5,0.3); //NB informative
  }
'


  model_stan <- paste(
    data_com,

    param_com,

    if(time_method[[1]] == 'ar') param_phi_ar,
    if(time_method[[1]] == 'var') param_phi_var,

    '}',

    #   tparam_com,
    if(method == 'lh')tparam_lh,
    if(method == 'eck')tparam_eck,

    model_com
  )

  `%notin%` <- Negate(`%in%`)


  # data part
  old <- options()        # old options
  on.exit(options(old)) 	# reset once exit the function

  options(na.action='na.pass') # to preserve the NAs

  ii <- rep(c(1:2),length.out = passes)

  out_list <- mylm(formula = formula, data = data) # produces the design matrix

  X <- out_list$X


  y_obs <- data[!is.na(data$y),]$y

  # index for observed values

  data$date_id <- as.numeric(factor(data$date))

  i_y_obs <- data[!is.na(data$y),]$date_id

  # index for missing values
  i_y_mis <- data[is.na(data$y),]$date_id



  data_list <- list(ii = ii, # 1 for odd/forward passes 2 for even/backward
                    passes = passes,
                    T = nrow(df), # time points
                    y_obs = y_obs,# y values in the obs df
                    N_y_obs = length(i_y_obs),  #nrow(i_y_obs) numb obs points
                    N_y_mis = length(i_y_mis), #nrow(i_y_mis) numb preds points

                    i_y_obs = i_y_obs, # index of obs points
                    i_y_mis = i_y_mis, # index of preds points
                    K = ncol(X),
                    X = X,
                    bfi = bfi
  ) # a list with all the distance/weights matrices

  ini <- function(){list(sigma = 1,
                         alpha =  0.9,
                         y =  rep(80,nrow(df)),
                         f =  matrix(rep(50, passes * nrow(df)), nrow = nrow(df),  ncol = passes),
                         y0 =  matrix(rep(70, passes * nrow(df)), nrow = nrow(df),  ncol = passes),
                         b =  matrix(rep(50, passes * nrow(df)), nrow = nrow(df),  ncol = passes)
  )}


  # print(model_stan)

  fit <- rstan::stan(model_code = model_stan,
                     model_name = "model_stan",
                     data = data_list,
                     #pars = pars,
                     iter = iter,
                     warmup = warmup,
                     init = ini,
                     chains = chains,
                     verbose = FALSE,
                     seed = seed,
                     refresh = refresh
  )
  attributes(fit)$formula <- formula

  #class(fit) <- 'ssnbayes'

  fit
}





#' A simple modeling function using a formula and data
#'
#' @param formula A formula as in lm()
#' @param data A data.frame containing the elements specified in the formula
#' @return A list of matrices
#' @importFrom stats model.matrix model.response
#' @export
#' @author Jay ver Hoef
#' @examples
#' options(na.action='na.pass')
#' data("iris")
#' out_list = mylm(formula = Petal.Length ~ Sepal.Length + Sepal.Width, data = iris)


mylm <- function(formula, data) {
  # get response as a vector
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data"), names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE
  mf[[1L]] <- as.name("model.frame")
  mf <- eval(mf, parent.frame())
  y <- as.vector(model.response(mf, "numeric"))
  # create design matrix
  X <- model.matrix(formula, data)
  # return a list of response vector and design matrix
  return(list(y = y, X = X))
}



