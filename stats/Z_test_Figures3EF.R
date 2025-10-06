# install.packages("brglm2")  # run once
library(brglm2)

# --- clonal level ---
# XDRTCCs
N1 <- 386 # total number of XDRTCCs defined
pre1 <- 30 # XDRTCCs clones detected PreTx
post1 <- 69 # XDRTCCs clones detected PostTx
# NonXDRTCCs
N2 <- 3635 # total number of NonXDRTCCs defined
pre2 <- 3635 # NonXDRTCCs clones detected PostTx 
post2 <- 480  # Non XDRTCCs clones detected PostTx

# Build aggregated binomial rows
df <- data.frame(
  detected = c(post1, pre1,  post2, pre2),
  notdet   = c(N1-post1, N1-pre1, N2-post2, N2-pre2),
  time     = factor(c("Post","Pre","Post","Pre")),
  group    = factor(c("XDR","XDR","Non","Non"))
)

# Set reference levels so the interaction is XDR vs Non, Post vs Pre
df$group <- relevel(df$group, ref = "Non")
df$time  <- relevel(df$time,  ref = "Pre")

# Firth-type bias-reduced logistic regression
fit <- glm(cbind(detected, notdet) ~ group * time,
           family = binomial(), data = df,
           method = "brglmFit")

# Extract the interaction (log Ratio-of-Odds-Ratios for XDR vs Non)
coefs <- coef(fit)
vc    <- vcov(fit)
beta  <- coefs["groupXDR:timePost"]
se    <- sqrt(vc["groupXDR:timePost","groupXDR:timePost"])

z_wald <- beta / se
p_two  <- 2 * pnorm(-abs(z_wald))
ROR    <- exp(beta)                           # Ratio of odds ratios (XDR vs Non)
CI95   <- exp(beta + c(-1.96, 1.96) * se)

cat(sprintf("Interaction (XDR vs Non, Post vs Pre):\n  Z = %.4f, p = %.4g\n  ROR = %.4f, 95%% CI = [%.4f, %.4f]\n",
            z_wald, p_two, ROR, CI95[1], CI95[2]))

##### Template level

# Counts
# same definitions as above, except pre1/2 and post1/2 are templates
N1 <- 386;  pre1 <- 109;  post1 <- 4423 
N2 <- 3635; pre2 <- 3920; post2 <- 5083

# If you have sequencing depth per stratum, put those in; else all 1
depth_pre_XDR  <- 1
depth_post_XDR <- 1
depth_pre_Non  <- 1
depth_post_Non <- 1

df <- data.frame(
  count = c(post1, pre1,  post2, pre2),
  group = factor(c("XDR","XDR","Non","Non"), levels = c("Non","XDR")),
  time  = factor(c("Post","Pre","Post","Pre"), levels = c("Pre","Post")),
  Ndef  = c(N1, N1, N2, N2),
  depth = c(depth_post_XDR, depth_pre_XDR, depth_post_Non, depth_pre_Non)
)

# Poisson log-linear model
fit <- glm(count ~ group * time + offset(log(Ndef)) + offset(log(depth)),
           family = poisson, data = df)

summary(fit)

# Extract interaction (XDR vs Non, Post vs Pre)
co <- summary(fit)$coefficients
beta <- co["groupXDR:timePost","Estimate"]
se   <- co["groupXDR:timePost","Std. Error"]
z    <- beta / se
pval <- 2 * pnorm(-abs(z))
RRR  <- exp(beta)  # ratio of rate ratios
CI95 <- exp(beta + c(-1.96, 1.96)*se)

cat(sprintf("Interaction (XDR vs Non, Post vs Pre):\n  Z = %.4f, p = %.4g\n  RRR = %.3f, 95%% CI = [%.3f, %.3f]\n",
            z, pval, RRR, CI95[1], CI95[2]))

