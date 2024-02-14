
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
	int N;     // number of surveys
	int Y;     // number of years
	int P;     // number of survey groups (partitions)
	int S[P];  // number of survey sites per partition
	
	// look-up vectors
	int XY[N];
	int XS[N];
	
	// survey count 
	// data
	int y[N]; 
	
	// survey area size
    vector[Y] survey_area[sum(S)];
}
transformed data {

	// total number 
    // of survey sites
	int sumS = sum(S);
		
    // construct
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
    
    // stability
    real<lower=0> gamma;
	
	// random effects error
    real<lower=0> tauS[P];
    real<lower=0> phi;
}
transformed parameters {

	// linear predictors
	vector[Y] xt_log[sumS];
	vector[Y] lambda[sumS];
    vector[Y] scalar_log[sumS];
	
	vector[sumS] alpha_re;
    
    vector[N] eta;
		
	{
		int k = 0;
		for(i in 1:P) {			
			for (j in 1:S[i]) {
			
				k += 1;
				
				// centred random effects
				alpha_re[k] = alphaP[i] + tauS[i] * alphaS[k];
				
				// define model for growth rate (lambda)
				lambda[k] = alpha_re[k] * year;
				
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
    
    // random effects
    alphaS ~ std_normal();
    
    // standard error terms
    square(tauS) ~ exponential(1);
	
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
    
	// POPULATION TRENDS
    {	
        int  k;	
        real tau_sum;
        
        // PER SITE
        k = 0;
		for (i in 1:P) {			
			for (j in 1:S[i]) {
			
				k += 1;
                trend_per_site[k]   = exp(alpha_re[k] * T);
                numbers_per_site[k] = exp(xt_log[k] + scalar_log[k]);
            }
        }
        
        // PER PARTITION
        for (i in 1:P) {		
        
            tau_sum = sqrt(T * square(tauS[i]));
            trend_per_partition[i] = exp(alphaP[i] * T + square(tau_sum) / 2.0);
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
	traces[3] = 0.0;
    traces[4] = vector_norm(alphaS);
	traces[5] = 0.0;
    traces[6] = vector_norm(to_vector(tauS));
    traces[7] = phi;
    traces[8] = vector_norm(trend_per_partition);
}
