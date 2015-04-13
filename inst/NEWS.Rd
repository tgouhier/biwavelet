\name{NEWS}
\title{News for Package 'biwavelet'}
\encoding{UTF-8}

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
