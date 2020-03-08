---
layout: default
mathjax: true
---

{% include mathjax.html %}

***

[Supplementary Methods](SM)

***

## Introduction 
Spread of COVID-19 in early 2020 has raised important concerns about the ability of national health systems to detect the positive cases, as well as the intervention rate a country is expected to put in place to contain the infection.

  Using data on desease spread and containment through the public repository [CSSEGISandData/COVID-19](https://github.com/CSSEGISandData/COVID-19) we established a predictive model to estimate the outbreak of the infection in a given population. 

In some regions of China the infection rate has significantly decreased compared to the initial exponential spread of the infection. This information can be used to build predictive models which can help other countries to estimate the extent of the outbreak.  
First we modeled the number of infections over time based on the desease outbreak in China. Our model captures the initial exponential phase of the outbreak and the effect of the external intervention to contain the infection. Thus the exponential phase termimnates with a peak of maximum rate of new infected cases and then is followed by a critical period where the number of infections can decrease if the intervention remains stable.  

|<img src="Figures/Figure_stat_1.png"/>|
|:--:|
|Fig. 1: Dynamics of the infection in the province of Henan, Mainland China|

## The model
We used a generative model to describe the dynamics of the infected population in a given geographic area. Our model takes into account the effect of the local interventions by coupling the average number of infections $$X(t)$$ with a dynamical variable $$A(t)$$ which acts against the spread. The observed number of cases is then generated from a Poisson distribution with rate $X(t)\cdot p(t)$, where $p(t)$ is a fraction of the true number of infections which varies over time (due for instance to the increased number of tests during the acute phase).

The model is defined by two differential equations for the size of the infected population $X(t)$ and the action $A(t)$ and a probabilistic rule on the observed infections 

$$
\begin{align*}
\frac{dX}{dt} &= \lambda X - X A,\quad X(t_0) = X_0\\
\frac{dA}{dt} &= h X(t)p(t),\quad A(t_0)=0\\
n(t) & \sim \mathrm{Poisson}(X(t) \cdot p(t)),\quad p(t)=\frac{X(t)^g}{X(t)^h+k^g}
\end{align*}
$$ 

where $\lambda$ is the (observed) infection rate and $h$ quantifies the effect of local interventions.  
$X_0$ is the average number of infections at the initial time (22/01/2020).
The lack of observations at early stages is captured by the factor $p(t)$ which depends directly on the daily infection rate. The assumption being that the more infections are detected the more tests and controls are put in place to monitor the infected population. Note however that this is not sufficient to estimate the true number of infected individual, but to accommodate the low number of observations when the epidemics reaches a given country.

Figure 2 illustrates the typical dynamics of the model.

|<img src="Figures/Figure_1.png"/>|
|:--:|
|Fig. 2: *Dynamics of the number of infections*|

## Statistical inference of model parameters
By using the available daily reported cases in the public repository [CSSEGISandData/COVID-19](https://github.com/CSSEGISandData/COVID-19) we can estimate the parameters of the model from the data for each country/region affected by the infection. Knowing the model parameters allow us to draw predictions on how the epidemics will evolve. For this analysis we assume that the rate of infection $\lambda$ is the same for all countries whereas all the other model parameters are country-dependent. This allows us to exploit the worldwide data to strengthen the predictive power of the model.

The framework of statistical inference allows us to estimate the model parameters and make predictions while taking into account statistical uncertainties derived from the data and the prior uncertainty. We performed a global analysis on 173 countries included in the CSSE dataset [1]. 

The interactive chart below gives an overview of the course of the infection for each country.
<iframe width="800" height="400" frameborder="0" scrolling="no"
src="notebook/plotly_chart.html"></iframe>

## References
1. Dong, Ensheng, Hongru Du, and Lauren Gardner. "An interactive web-based dashboard to track COVID-19 in real time." The Lancet Infectious Diseases (2020).   

