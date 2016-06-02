\name{NEWS}
\title{News for Package 'biwavelet'}
\encoding{UTF-8}

\section{Changes in biwavelet version 0.20.7 (2016-06-01)}{
  \subsection{fixed}{
    \itemize{
      \item Fixed handling of time in \code{check.data}.
    }
  }
}

\section{Changes in biwavelet version 0.20.6 (2016-06-01)}{
  \subsection{fixed}{
    \itemize{
      \item Fixed x-axis in \code{plot.biwavelet}.
    }
  }
}

\section{Changes in biwavelet version 0.20.3 (2016-05-08)}{
  \subsection{fixed}{
    \itemize{
      \item Fixed displacement of COI, contours and phase arrows in \code{plot.biwavelet} when adding a color bar.
    }
  }
}

\section{Changes in biwavelet version 0.20.2 (2016-05-06)}{
  \subsection{fixed}{
    \itemize{
    \item Fixed \code{check.datum}; hopefully for the last time.
    \item Faster wt bases and row quantile (Rcpp implementations):
      \itemize{
      \item The \code{param} parameter for all \code{rcpp_wt_bases_*} must be
            within interval (0..10).
      \item The \code{rcpp_row_quantile} function requires
            a matrix as a parameter (use \code{as.matrix()} for vectors).
      }
    }
  }
}

\section{Changes in biwavelet version 0.19.2 (2016-05-03)}{
  \subsection{fixed}{
    \itemize{
    \item Fixed Rcpp implementation of the \code{wt.bases} functions,
    i.e., \code{rcpp_wt_bases_*.cpp}.
    Replacing \code{int} type with \code{double} type for \code{scale} parameter
    which caused small scales to be rendered incorrectly.
    }
  }
}
\section{Changes in biwavelet version 0.19.1 (2016-04-29)}{
  \subsection{fixed}{
    \itemize{
    \item Fixed interpretation of phase differences in \code{plot.biwavelet} help file
    }
  }
  \subsection{changed}{
    \itemize{
      \item Added unit tests for approx 78\% of the code.
      \item Implemented a parallelized Monte Carlo simulation function
            \code{wtc_sig_parallel} which is 2 to 4 times faster on a 4-core CPU
            than the original \code{wtc.sig}. The speedup is noticeable on:
            \enumerate{
              \item large simulations \code{nrads >= 800},
              \item multiple simulations,
              \item multi-core systems with 4+ cores.
            }
        However, parallelization involves a significant heat-up phase because all
        the workers need to be started and they need to load all the required
        packages. This will be addresses in future versions of biwavelet.
      \item Added a speed-optimized version of \code{convolve2D}.
      \item Replaced standard \code{arima.sim} function with a pair of functions
        \code{get_minroots} and \code{ar1_ma0_sim}. These functions are still
        implemented in R. We can reimplement them later in C.
      \item Reimplemented \code{wt.bases} morlet, paul and dog in C.
      \item Removed unused function \code{meshgrid}.
    }
  }
}

\section{Changes in biwavelet version 0.17.11 (2015-10-09)}{
  \subsection{fixed}{
    \itemize{
    \item close all progress bars after use
    \item Function \code{wtc} can now handle non-finite values when computing
    the quantiles of the rsq values from the Monte Carlo simulations
    }
  }
}

\section{Changes in biwavelet version 0.17.10 (2015-04-29)}{
  \subsection{fixed}{
    \itemize{
    \item Added ability to handle custom color palettes in \code{plot.biwavelet}.
    Users can now specify any color scheme using the \code{fill.cols} argument.
    }
  }
}

\section{Changes in biwavelet version 0.17.9 (2015-04-29)}{
  \subsection{fixed}{
    \itemize{
    \item Fixed limited padding issue, which could lead to weird edge effects.
    Current padding level is identical to that of Torrence & Compo (1998).
    \item Changed the default \code{tol} value from 0.95 to 1 in the \code{plot} function.}
    }
}

\section{Changes in biwavelet version 0.17.8 (2015-04-28)}{
  \subsection{fixed}{
    \itemize{
    \item Added semi-transparent COI.}
    }
}

\section{Changes in biwavelet version 0.17.7 (2015-04-13)}{
  \subsection{fixed}{
    \itemize{
    \item Fixed \code{check.datum} function so that it does not assume a sampling frequency of 1.}
    }
}

\section{Changes in biwavelet version 0.17.6 (2015-04-05)}{
  \subsection{fixed}{
    \itemize{
    \item Added ability to set \code{zlim} in \code{plot.biwavelet}.}
    }
}

\section{Changes in biwavelet version 0.17.5 (2014-11-05)}{
  \subsection{fixed}{
    \itemize{
    \item Improved the implementation of \code{phase.plot} to allow for much better looking phase arrows (thanks Huidong Tang).}
    }
}

\section{Changes in biwavelet version 0.17.4 (2014-11-04)}{
  \subsection{fixed}{
    \itemize{
    \item Made function \code{wt} faster by avoiding excessive padding (thanks Huidong Tang).}
    }
}

\section{Changes in biwavelet version 0.17.3 (2014-04-27)}{
  \subsection{fixed}{
    \itemize{
    \item Made \code{check.datum} tolerate slight inconsistencies in the size of timesteps.}
    }
}

\section{Changes in biwavelet version 0.17.2 (2014-04-11)}{
  \subsection{fixed}{
    \itemize{
    \item Added arguments in \code{plot.biwavelet} and \code{phase.plot} to control the length
      of the phase arrows and the size of the arrow heads independently.}
    }
}

\section{Changes in biwavelet version 0.17.1 (2013-07-15)}{
  \subsection{fixed}{
    \itemize{
    \item Fixed code in \code{check.data} to test for constant step size in the data.}
    }
}

\section{Changes in biwavelet version 0.17 (2013-06-05)}{
  \subsection{added}{
    \itemize{
      \item Function \code{pwtc} can be used to perform partial wavelet coherence between two time series
      \code{y} and \code{x1} by controlling for (or partialling-out) a third time series \code{x2}.
    }
  }
}
\section{Changes in biwavelet version 0.16 (2013-05-07)}{
  \subsection{added}{
    \itemize{
      \item Users can now specify the density of the phase arrows using the \code{plot} function.
    }
  }
  \subsection{fixed}{
    \itemize{
      \item Fixed bug in \code{wt} affecting the significance region (thanks Patrick Kilduff and Flora Cordoleani).
    }
  }
}
\section{Changes in biwavelet version 0.15 (2013-04-08)}{
  \subsection{added}{
    \itemize{
      \item Users can now specify the color, line width and line type for
      the COI, significance contours and phase arrows using the \code{plot} function.
    }
  }
    \subsection{fixed}{
\itemize{
      \item Removed misleading examples showing how to compute the 'bias-corrected' wavelet coherence.
      There is no bias for the wavelet coherence function, so using the default \code{type} argument in the
      \code{plot} function is recommended.
      \item Fixed typos in the documentation of plot.biwavelet and xwt (thanks Lei Cheng).
    }
  }
}
\section{Changes in biwavelet version 0.14 (2013-03-06)}{
  \subsection{added}{
    \itemize{
      \item As of biwavelet version 0.14, the bias-corrected wavelet and
        cross-wavelet spectra are automatically computed and plotted by default
        using the methods described by Liu et al. (2007) and Veleda et al. (2012).
        This correction is needed because the traditional approach for computing the power
        spectrum (e.g., Torrence and Compo 1998) leads to an artificial and systematic
        reduction in power at lower periods.
    }
  }
  \subsection{fixed}{
    \itemize{
      \item Plotting function now accepts traditional plotting flags such as xaxt and yaxt
        to control x and y tickmarks.
    }
  }
}
