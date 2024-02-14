
// NOTATION
// DIMENSIONS
// P: number of partitions
// N: number of records
// S: number of sites per partition
// Y: number of years
// LOOK-UP VECTORS
// XY: year
// XS: survey site
// OBSERVATIONAL DATA
// y: survey numbers

functions {

	real vector_norm(vector x) {
	    
	    real i = 0.0;
	    
	    for (j in 1:num_elements(x))
	        i += pow(x[j], 2.0);
	        
	    return pow(i, 0.5);
	}
}
data {
	
	// dimensions
	int N;
	int Y;
	int P; 
	int S[P];
	
	// look-up vectors
	int XY[N];
	int XS[N];
	
	// site level 
	// survey count 
	// data
	int y[N]; 
	
	// covariates
    vector[Y] survey_area[sum(S)];
}
transformed data {

	// number of sites
	int sumS = sum(S);
		    
    // year vector        
	vector[Y] year;
	for(i in 1:Y) {
		year[i] = i - 1;
	}
}
parameters {
	
    // density intercept parameters
	vector[sumS] x0_log;
	
	// mortality parameters
    vector[P]    alphaP;
	vector[sumS] alphaS;
    
    vector[P]    betaP;
	vector[sumS] betaS;
    
    // stability
    real<lower=0> gamma;
	
	// random effects error
    real<lower=0> tauS[2, P];
    real<lower=0> phi;
}
transformed parameters {
	
	// linear predictors
	vector[Y] xt_log[sumS];
	vector[Y] lambda[sumS];
    vector[Y] scalar_log[sumS];
	
	vector[sumS] alpha_re;
    vector[sumS] beta_re;
    
    vector[N] eta;
		
	{
		int k = 0;
		for(i in 1:P) {			
			for (j in 1:S[i]) {
			
				k += 1;
				
				// centred random effects
				alpha_re[k]  = alphaP[i] + tauS[1, i] * alphaS[k];
                beta_re[k]   = betaP[i]  + tauS[2, i] * betaS[k];
				
				// define model for growth rate (lambda)
				lambda[k] = (alpha_re[k] + beta_re[k] * year) .* year;
				
				// define model for density
				xt_log[k] = x0_log[k] + lambda[k];
                
                // define area scalar
                scalar_log[k] = gamma * log(survey_area[k]);
			}
		}
	}
    
    for (i in 1:N) {
        eta[i] = xt_log[XS[i]][XY[i]] + scalar_log[XS[i]][XY[i]];
    }
}
model {
		
	// likelihood
	y ~ neg_binomial_2_log(eta, phi);
	    
    // intercept values
    alphaP ~ normal(0, 1000);
    betaP  ~ normal(0, 1000);
    
    // random effects
    alphaS ~ std_normal();
    betaS  ~ std_normal();
    
    // standard error terms
    square(tauS[1]) ~ exponential(1);
    square(tauS[2]) ~ exponential(1);
    
    // stability
    gamma ~ gamma(1.0, 1.0);
    
    // inverse-dispersion
    phi ~ std_normal();
}
generated quantities {
	
	real traces[8];
    
    vector[Y]    numbers_per_site[sumS];
    vector[sumS] trend_per_site;
    vector[P]    trend_per_partition;
    
    vector[N] y_hat;
    vector[N] y_sim;
	vector[N] log_lik;
    
    // assessment time frame
    real T = Y - 1;
    real Tplus = T * (T + 1) / 2.0;
    
	// POPULATION TRENDS
    {	
        int  k;	
        real tau_sum;
        
        // PER SITE
        k = 0;
		for (i in 1:P) {			
			for (j in 1:S[i]) {
			
				k += 1;
                trend_per_site[k]   = exp(alpha_re[k] * T + beta_re[k] * square(T));
                numbers_per_site[k] = exp(xt_log[k] + scalar_log[k]);
            }
        }
        
        // PER PARTITION
        for (i in 1:P) {		
            
            // variance for linear term
            tau_sum = T * square(tauS[1, i]);
            // plus variance for second order term
            for (j in 0:(Y - 1)) {
                tau_sum += square(j * tauS[2, i]);
            }
            // standard error of the sum
            tau_sum = sqrt(tau_sum);
            // partition decline
            trend_per_partition[i] = exp(alphaP[i] * T + betaP[i] * Tplus + square(tau_sum) / 2.0);
        }
    }
    
    for (i in 1:N) {
    
        // expectation of the survey numbers data
		y_hat[i] = exp(eta[i]);
        
        // posterior prediction of the survey numbers data
		y_sim[i] = neg_binomial_2_log_rng(eta[i], phi);
        
        // log_likelihood for model comparison
        log_lik[i] = neg_binomial_2_log_lpmf(y[i] | eta[i], phi);
	}
    
	// SUMMARY STATISTICS FOR TRACE DIAGNOSTICS
	traces[1] = vector_norm(x0_log);
	traces[2] = vector_norm(alphaP);
	traces[3] = vector_norm(betaP);
    traces[4] = vector_norm(alphaS);
	traces[5] = vector_norm(betaS);
    traces[6] = vector_norm(to_vector(tauS[1])) + vector_norm(to_vector(tauS[2]));
    traces[7] = phi;
    traces[8] = vector_norm(trend_per_partition);
}
