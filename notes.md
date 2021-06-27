## COMPAS

### Misc Info

Florida state,  2013 - 2014, 18.6k people 

Most data from pretrial, was used to determine whether to release or not

On average, defendants in our dataset were not incarcerated for 622.87 days (sd: 329.19). ?



### Score

3 different risk scores: recidivism, violence and failure to appear (only use first 2)

“scores in the medium and high range garner more interest from supervision agencies than low scores, 
as a low score would suggest there is little risk of general recidivism,” 
Low: 1 - 4
Med: 5 - 7
High: 8 - 10



### Bias

Our logistic model found that the most predictive factor of a higher risk score was age. 
Defendants younger than 25 years old were 2.5 times as likely to get a higher score than middle aged offenders

While Black defendants had higher recidivism rates overall,  when adjusted for this difference and other factors, they were 45 percent more likely to get a higher score than whites.

Surprisingly, given their lower levels of criminality overall, female defendants were 19.4 percent more likely to get a higher score than men, controlling for the same factors.