---
layout: default
title: "Supplementary methods"
mathjax: true
permalink: /SM/
---

{% include mathjax.html %}

# Supplementary methods

## Data
The public repository [CSSEGISandData/COVID-19](https://github.com/CSSEGISandData/COVID-19) contains the list of confirmed, recovered and death cases updated every day. To get the number of infected individual at each day we consider the difference between confirmed and recovered or deats cases.

## Likelihood calculation
The data likelihood is calculated as

$$
L = \prod_{i\in countries}\prod_{t\in times}\mathrm{Pois}(I_{it}|X_{it}(h,x_0,k,g))
$$

where $X_{it}(h,x_0,k,g)$ is the average number of infected cases in country $i$ at time $t$, which was obtained by numerical integration of the differential equations of the model. 

## Prior distributions
As prior distributions we used gamma distributions as follows

$$
\begin{align}
\lambda&\sim\mathrm{Gamma}(1,10)\\
h&\sim\mathrm{Gamma}(5,100)\\
x_0&\sim\mathrm{Gamma}(1,0.1)\\
k&\sim\mathrm{Gamma}(4,0.01)\\
g&\sim\mathrm{Gamma}(1,1)
\end{align}
$$

## Sampling strategy
To draw samples from the posterior distribution of the model parameters we employed a standard Metropolis-Hastings scheme with a proposal obtained by combining a gamma distribution with equal scale and shape (mean at 1 to allow for local sampling around the current value) and the corresponding prior to allow for larger transitions:  

$$ Q(\theta^*,\theta)=\frac{1}{2}\mathrm{Gamma}(\theta^*/\theta,a,b)+\frac{1}{2}\pi(\theta)
$$

## Supplementary figures
<img src="../Figures/Figure_stat_2.png"/>
<imgcaption>Supplementary figure 1: parameter $h$ per country</imgcaption>

<img src="../Figures/Figure_stat_3.png"/>
<imgcaption>Suppementary figure 2: infection rate the 22nd of Jan 2020</imgcaption>

<img src="../Figures/Figure_stat_4.png"/>
<imgcaption>Suppementary figure 3: Hill scale parameter $k$ per country. This parameter can be interpreted as the size of the infected population at which the intervention becomes more pronounced.</imgcaption>

<img src="../Figures/Figure_stat_5.png"/>
<imgcaption>Suppementary figure 4: Hill shape parameter $g$ per country. Combined with the Hill scale, this parameter characterizes the strength of the intervention. In particular, after the infected population has surpassed the scale $k$, higher values of the shape $g$ denote a more pronounced switch in the actions taken to contain the spread.</imgcaption>

<img src="../Figures/Figure_stat_lambda.png"/>
<imgcaption>Suppementary figure 5: Posterior distribution of the daily infection rate $\lambda$</imgcaption>


