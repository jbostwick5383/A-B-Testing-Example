
install.packages("tidyverse")

#loading packages
library(tidyverse)

#reading in the dataset
AB_Test <- readxl::read_excel("website_data.xlsx")

View(AB_Test)

#Filter out conversions for variant_A
variant_a <- AB_Test %>% 
  filter(variant == "A" & converted == "TRUE")

#Total number of conversions for variant_A
conversions_A <- nrow(variant_a)

#Total number of visitors for variant_A
visitors_a <- nrow(AB_Test %>% filter(variant == "A"))

#conversion rate for variant_a
conversion_rate_a <- conversions_A/visitors_a

####Conversion rate B####

variant_b <- AB_Test %>%
  filter(variant == "B" & converted =="TRUE")

conversions_B <- nrow(variant_b)

#Total number of visitors for variant_b
visitors_b <- nrow(AB_Test %>% filter(variant == "B"))

#conversion rate for variant_b
conversion_rate_b <- conversions_B/visitors_b

#conversion rates
conversion_rate_a
conversion_rate_b

#computing relative uplift of conversion rates A & b
#the uplift is the % increase
uplift <- (conversion_rate_b - conversion_rate_a)/ conversion_rate_a * 100
uplift
#B is better than A by 83%

####Computing the pooled probability, standard error, the margin of error, and 
#difference in proportion (point estimate) for variants A & B

#pooled sample proportion
p_pool <- (conversions_A + conversions_B)/(visitors_a + visitors_b)
print(p_pool)

#compute the standard error for variants A & B
se_pool <- sqrt(p_pool * (1-p_pool)*((1/visitors_a)+(1/visitors_b)))

print(se_pool)

#compute the margin of error for the pool
MOE <- se_pool * qnorm(0.975)

print(MOE)

#point estimate or difference in proportion
d_hat <- conversion_rate_b - conversion_rate_a
print(d_hat)

####computing the z score####
#compute the z score so that we can determine the p-value
z_score <- d_hat/se_pool
print(z_score)

#now we can use the z score to determine the p-value
p_value <- pnorm(q = -z_score, mean = 0, sd = 1) *2
print(p_value)

#now compute the confidence interval for the pool
#lets compute confidence interval for the pool using pre-calculated results
ci <- c(d_hat - MOE, d_hat + MOE)
print(ci)

#Using the same steps now we will compute the confidence interval for variants A seperately
x_hat_a <- conversions_A/visitors_a
se_hat_a <- sqrt(x_hat_a * (1- x_hat_a)/visitors_a)
ci_a <- c(x_hat_a - qnorm(0.975) * se_hat_a, x_hat_a + qnorm(0.975) * se_hat_a)
print(ci_a)

#using the same steps as already shown, conput the confidence interval for variants b seperately
x_hat_b <- conversions_B/visitors_b
se_hat_b <- sqrt(x_hat_b * (1 - x_hat_b)/visitors_b)
ci_b <- c(x_hat_b - qnorm(0.975) * se_hat_b,
          x_hat_b + qnorm(0.975) * se_hat_b)
print(ci_b)

#visualize the results computed so far in a dataframe
vis_result_pool <- data.frame(
  metric = c(
    'Estimated Difference',
    'Relative Uplift(%)',
    'pooled sample proportion',
    'Standard Error of Difference',
    'z-score',
    'p-value',
    'Margin of Error'),
  value = c(
    conversion_rate_b - conversion_rate_a,
    uplift,
    p_pool,
    se_pool,
    z_score,
    p_value,
    MOE)
)

#Recommendations & Conclusions
# - Variant A has 20 conversions and 721 hits whereas Variant B has 37 conversions and 730 hits.
# - Relative uplift of 82.72% based on a variant A conversion rate is 2.77% and for B is 5.07%.  Hence, variant B is better than A by 82.72%.
# - For this analysis P-value computed was 0.02448.  Hence, there is strong statistical significance in the test results.
# - From the above results that depict strong statistical significance.  You should reject the null hypothesis and proceed with the launch.
# - Therefore, Accept Variant B and you can roll it to the users for 100%.

#Limitations
#It is one of the tools for conversion optimization and its not an independent solution 
#and its not an independent solution and its not going to fix all the conversion issues
#of ours and it cant fix the issues as you get with messy data and you need to perform more
#than just an A/B test to improve on conversions.
