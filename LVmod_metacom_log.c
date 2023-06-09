/* File LVmod_metacom_log.c */
#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
static double *Pars;

/* Initializer  */
void initmod(void(* odeparms)(int *, double *)) {
  DL_FUNC get_deSolve_gparms;
  SEXP gparms;
  get_deSolve_gparms = R_GetCCallable("deSolve","get_deSolve_gparms");
  gparms = get_deSolve_gparms();
  Pars = REAL(gparms);
}

/* Derivatives */
void derivs (int *neq, double *t, double *y, double *ydot,
             double *yout, int *ip) {
  const int cadj = 0;  // Indexing parameter. Allows easier transfer from C (indexes from 0) to R (indexes from 1).
  int N = Pars[0+cadj];
  int M = Pars[1+cadj];
  int i;
  int m;
  int minval = Pars[2+cadj];
  double dl = 1/Pars[3+N*(2+N)+cadj+4]; // disturbance rate
  double p_incr = 1/(double)(M);

  // disturbance function
  ydot[N*M+cadj] = -y[N*M+cadj]*dl; // waiting time function
  ydot[N*M+1+cadj] = 0; // dummy variable, storing random uniform number

  m = 0;
  while(m < M) {
    // all species in patch m
    i = 0;
    while(i < N) {
      if(m==0) {
        yout[i] = 0; // re-set
      }
      if(y[i+N*m+cadj]<=minval) {
      	ydot[i+N*m+cadj] = 0;
      } else {
        double comp_force = 0;
        int j = 0;
        while(j < N) {
          double tmp = y[j+N*m+cadj];
          if(tmp <= minval) {
            tmp = 0;          
          } else {
            tmp = exp(tmp);
          }
        
          comp_force = comp_force + tmp*Pars[3+cadj+2*N+i+N*j];
          
          j = j+1;
        }
        ydot[i+N*m+cadj] = Pars[3+i+cadj+N]*(Pars[3+i+cadj]+comp_force);
        yout[i] = yout[i]+p_incr; // p
      }
      i = i+1;
    }
    m = m+1;
  }

  // dispersal rates
  i = 0;
  while(i < N) {
    // rate is c*p*(1-p)
    double cpp = Pars[3+N*(2+N)+cadj+5+i]*yout[i]*(1-yout[i]);

    ydot[N*M+2+i+cadj] = -y[N*M+2+i+cadj]*cpp; // waiting time function
    ydot[N*M+2+N+i+cadj] = 0; // dummy variable, storing random uniform number

    i = i+1;
  }

}

/* Events */
void event(int *n, double *t, double *y) {
  GetRNGstate();
  double ysd = 0;
  const int cadj = 0;
  int N = Pars[cadj];
  int M = Pars[1+cadj];
  int minval = Pars[2+cadj];
  int i = 0;
  int m = 0;
  double tmp = 0;

  double dc = Pars[3+N*(2+N)+cadj];
  double dz = Pars[3+N*(2+N)+cadj+1];
  double dn = Pars[3+N*(2+N)+cadj+2];
  double dm = Pars[3+N*(2+N)+cadj+3];

  // find next event
  int which_event = 0;
  double min_event = 1;
  i = 0;
  while(i < (N+1)) {
    if(i==0) {
      tmp = (y[N*M+1+cadj]-(1-y[N*M+cadj]));
    } else {
      tmp = y[N*M+2+N+(i-1)+cadj]-(1-y[N*M+2+(i-1)+cadj]);
    }
    if(tmp < min_event) {
      which_event = i;
      min_event = tmp;
    }
    i = i+1;
  }

  // test disturbance trigger
  if(which_event == 0) {
    // if disturbance...
    // re-set the disturbance counter
    y[N*M+1+cadj] = unif_rand(); // new uniform function
    y[N*M+cadj] = 1; // reset to initial value

    // implement disturbance
    m = 0;
    while(m < M) {
      i = 0;
      while(i < N) {
        //R_CheckUserInterrupt();
        if(y[i+N*m+cadj] >= minval) {
          ysd = sqrt(dc*pow(exp(y[i+N*m+cadj]),dz));
          double tmp = exp(y[i+N*m+cadj]) + (ysd+dn)*norm_rand()+dm;
          if(tmp > 0) {
            y[i+N*m+cadj] = log(tmp);
          } else {
            y[i+N*m+cadj] = minval;
          }
        }
        i = i+1;
      }
      m = m+1;
    }
  } else {
    // if dispersal...
    i = which_event-1; // species i is the dispersing species

    // find and choose empty sites to disperse into
    int startpos = floor(unif_rand()*M);
    int dispos = 0;
    int m = 0;
    while(m < M) {
      dispos = (m+startpos)%M+cadj;
      if(y[i+N*dispos+cadj]<=minval) { // site is empty
        m = M; // end loop
      } else {
        m = m+1;
      }
    }  // species is dispersing into site dispos

    // increase y for dispersed species and site
    y[i+N*dispos+cadj] = log(Pars[3+N*(3+N)+cadj+5+i]);

    // re-set the dispersal counter
    y[N*M+2+N+i+cadj] = unif_rand(); // new uniform function
    y[N*M+2+i+cadj] = 1; // reset to initial value
  }
  PutRNGstate();
}

void myroot(int *neq, double *t, double *y, int *ng, double *gout, double *out, int *ip) {
  const int cadj = 0;
  int N = Pars[0+cadj];
  int M = Pars[1+cadj];

  // check for disturbance
  double tmp = y[N*M+1+cadj]-(1-y[N*M+cadj]);
  //(tmp<0) {
  //  tmp = 0;
  //}
  gout[0] = tmp;

  // if not disturbance, then check for dispersal
  int i = 0;
  while(i < N) {
    double tmp = y[N*M+2+N+i+cadj]-(1-y[N*M+2+i+cadj]);
    //if(tmp<0) {
    //  tmp = 0;
    //}
    gout[i+1] = tmp;
    i = i+1;
  }
}

/* The Jacobian matrix */
/*
void jac(int *neq, double *t, double *y, int *ml, int *mu,
  double *pd, int *nrowpd, double *yout, int *ip) {
  pd[0] = -k1;
  pd[1] = k1;
  pd[2] = 0.0;
  pd[(*nrowpd)] = k2*y[2];
  pd[(*nrowpd) + 1] = -k2*y[2] - 2*k3*y[1];
  pd[(*nrowpd) + 2] = 2*k3*y[1];
  pd[(*nrowpd)*2] = k2*y[1];
  pd[2*(*nrowpd) + 1] = -k2 * y[1];
  pd[2*(*nrowpd) + 2] = 0.0;
}
*/

/*
void jac(int *neq, double *t, double *y, int *ml, int *mu,
  double *pd, int *nrowpd, double *yout, int *ip) {
  int cadj = 0;
  int N = Pars[0+cadj];
  
  int i = 0;
  while(i<N) {
    int j = 0;
    while(j<N) {
      if(i != j) {
        pd[i+j*N+cadj] = Pars[i+1+N+cadj]*y[i+cadj]*Pars[1+2*N+i+N*j+cadj];
      } else {
        pd[i+j*N+cadj] = Pars[i+1+N+cadj]*Pars[i+1+cadj]+
          2*Pars[i+1+N+cadj]*y[i+cadj]*Pars[1+2*N+i+N*j+cadj];
        int k = 0;
        while(k<N) {
          if(i!=k) {
            pd[i+j*N+cadj] = pd[i+j*N+cadj]+Pars[i+1+N+cadj]*y[k+cadj]*Pars[1+2*N+i+N*k+cadj];
          }
          k = k+1;
        }
      }
      j = j+1;
    }
    i = i+1;
  }
}
*/

/* END file LVmod_SEXP.c */
