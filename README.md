# Overview

This code is a test example for the multi-model forecast Kalman filter, as initialliy formulated in November 2020. In particular, we assume linear, Gaussian models and utilize some linear combination of our models to create a psuedo-truth signal for each model that that model observes. The purpose of this code, at least initially, is to be able to run the Gershgorin-Majda 2010 system, as well as the additive and multiplicative models.

# The Gershgorin-Majda 2010 System

Although this is not the official name for this test system, it was first proposed by Gershgorin and Majda in their 2010 paper *Test models for improving filtering with model errors through stochastic parameter estimation*. Specifically, the exactly solvable test model is

$$\begin{aligned}
	\frac{d\,u\left( t \right)}{dt} &= \left( -\gamma\left( t \right) + i\,\omega \right)\,u\left( t \right) + b\left( t \right) + f\left( t \right) + \sigma\,\dot{W}\left( t \right), \\
	\frac{d\,b\left( t \right)}{dt} &= \left( -\gamma_{b} + i\,\omega_{b} \right)\,\left( b\left( t \right) + \widehat{b} \right) + \sigma_{b}\,\dot{W}_{b}\left( t \right),\\
	\frac{d\,\gamma\left( t \right)}{dt} &= -d_{\gamma}\,\left( \gamma\left( t \right) + \widehat{\gamma} \right) + \sigma_{\gamma}\,\dot{W}_{\gamma}\left( t \right),
\end{aligned}$$

where $$u\left( t \right)$$ and $$b\left( t \right)$$ are complex-valued and $$\gamma\left( t \right)$$ is real-valued. The terms $$b\left( t \right)$$ and $$\gamma\left( t \right)$$ represent additive and multiplicative bias corrections terms. Also, $$\omega$ is the oscillation frequency of $$u\left( t \right)$$, $$f\left( t \right)$$ is external forcing, and $$\sigma$ characterizes the strength of the white noise forcing $$\dot{W}\left( t \right)$$. The parameters $$\gamma_{b}$$ and $$d_{\gamma}$$ represent the damping and parameters $$\sigma_{b}$$ and $$\sigma_{\gamma}$$ represent the strength of the white noise forcing of the additive and multiplicative bias correction terms, respectively. The parameters $$\widehat{b}$$ and $$\widehat{\gamma}$$ are the stationary mean bias correction values of $$b\left( t \right)$$ and $$\gamma\left( t \right)$$, correspondingly, and $$\omega_{b}$$ is the frequency of the additive noise. Note that the white noise $$\dot{W}_{\gamma}\left( t \right)$$ is real-valued while the white noises $$\dot{W}\left( t \right)$$ and $$\dot{W}_{b}\left( t \right)$$ are complex-valueds and their real and imaginary parts are independent real-valued white noises.

In Gershgorin and Majda's original paper, they consider an oscillatory external forcing

$$\begin{aligned}
	f\left( t \right) &= A_{f}\,e^{i\,\omega_{f}\,t}.
\end{aligned}$$

For the purpose of testing some parameter estimation techniques, Gershgorin and Majda considered an additive model of the original Gershgorin-Majda 2010 system where there is only additive bias correction

$$\begin{aligned}
	\frac{d\,u\left( t \right)}{dt} &= \left( -\overline{d} + i\,\omega \right)\,u\left( t \right) + b\left( t \right) + f\left( t \right) + \sigma\,\dot{W}\left( t \right), \\
	\frac{d\,b\left( t \right)}{dt} &= \left( -\gamma_{b} + i\,\omega_{b} \right)\,\left( b\left( t \right) + \widehat{b} \right) + \sigma_{b}\,\dot{W}_{b}\left( t \right),
\end{aligned}$$

where $$\overline{d}$$ is the mean value of the damping, as well as a multiplicative model where there is only multiplicative bias correction

$$\begin{aligned}
	\frac{d\,u\left( t \right)}{dt} &= \left( -\gamma\left( t \right) + i\,\omega \right)\,u\left( t \right) + f\left( t \right) + \sigma\,\dot{W}\left( t \right), \\
	\frac{d\,\gamma\left( t \right)}{dt} &= -d_{\gamma}\,\left( \gamma\left( t \right) + \widehat{\gamma} \right) + \sigma_{\gamma}\,\dot{W}_{\gamma}\left( t \right).
\end{aligned}$$
