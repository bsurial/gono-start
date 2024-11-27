# Sample Size calculations for my GonoStart trial

## Learnings: 

**26.11.2024 TIL from Andi:**  
When you want to power a non-inferiority study, you can use estimators for superiority trials and approximate the sample size by doing the following: 
Let's say we want a non-inferiority trial with a `delta of 1.5` (on a relative scale), `alpha of 0.025` (one-sided), and `power of 0.8`. You then take the superiority calculator and do the following: 

- effect-size = `delta (1.5)`
- one-sided alpha = `1 - Power (1 - 0.8 = 0.2)`
- Power = `1 - one-sided alpha (1 - 0.025 = 0.975)`

This should give you an approximation. I used simulation in `01-samplesized_simulation.R` but could have also used the method above. 


## Potentially eligible patients: 
![selection_of_sample](https://github.com/user-attachments/assets/94b34fa4-f963-4dff-84b6-2ba908fb32d7)
