#include <iostream>      // pour les sorties et entrées
#include <cmath>         // Fonctions mathématiques
#include <cstdlib>		 // Pour utiliser rand()
#include <vector>

using namespace std;

using Vector = vector<double>;


float gamma(int n, float C, float beta)
{
	return C/(pow(n+1,beta));
}

float F1(float x, float g){
	float p=0.95;
	if(g<=0){
		return(1-p);
	}
	else{
		return(-1*p);
	}
}
float F2(float x, float g)
{
    return x-g;
}

float gener_gauss()
{
	float U = rand()/(double)RAND_MAX ;
	float V = rand()/(double)RAND_MAX ;
	
	return sqrt(-2*log(U))*cos(2*M_PI*V);
}


void alogrithme_rm(float x0, int n, float C, float beta,float (*F)(float, float))
{
	float X[n];
	X[0] = x0;
	float G[n];
	
	for (int i = 0; i < n; i++)
	{
		G[i] = gener_gauss();
	}
	for (int i = 0; i < n-1; i++)
	{
		X[i+1] = X[i] - gamma(i,C,beta)*F(X[i],G[i+1]);
	}
	
	for (int i = 0; i < n-1; i++)
	{
		cout << X[i] << endl;
	}
}




void algorithme_rm2(float x0, int n, float C, float beta,float (*F)(float, float))
{
    double X[n];
	double G[n];
	X[0] = x0;
	
	for (int i = 0; i < n; i++)
	{
		G[i] = gener_gauss();
	}
	
    for (int i = 0; i < n-1; i++)
	{
        X[i+1] = X[i] - pow(gamma(i,C,beta),2)*F(X[i],G[i+1]);
	}
    for (int i = 0; i < n-1; i++)
	{
		cout << X[i] << endl;
	}
}
	

float F3(float x, float g)
{
    if (g<=x) return 1-0.95;
	return -0.95;
}


float alpha_star(float x_star)
{
	float c_star = exp(-pow(x_star,2)/2)/sqrt(2*3.14);
	return 1/(2*c_star);
}


double mean(Vector V)
{
	double res = 0;
	
	for (int i = 0; i < V.size(); i++)
	{
		res = res + V[i];
	}
	return res/V.size() ;
}

double max_element(double a, double b)
{
	if (a>=b) return a;
	return b;
}

double max_vect(Vector G)
{
	double m = G[0];
	
	for (int i = 0; i < G.size(); i++)
	{
		if (G[i] >= m) m = G[i];
	}
	
	return m;
}

double f(int m, Vector g, float r, float sigma, float T, float x, float K)
{
	int d = g.size();
    	float S_T[d];
	float res[d];
	
	for (int i = 0; i < d; i++) 
		S_T[i] = x*exp((r - sigma*sigma/2)*T + sigma*sqrt(T)*g[i]);
	for (int i = 0; i < d; i++)
	{
		res[i] = exp(-r*T)*max_element(S_T[i] - K, 0);
	}
	return res[m];
}


void vect_gauss(int n)
{
	double G[n];
	
	for (int i = 0; i < n; i++) G[i] = gener_gauss();
}


double approx_variance(double Lambda, float r, float sigma, float T, float x, float K, Vector G)
{
	Vector vect1;
	Vector vect2;
	
	for (int i = 0; i < G.size(); i++) 
	{
		vect1.push_back(exp(-Lambda*G[i]+0.5*Lambda*Lambda)*pow(f(i, G, r, sigma, T, x, K), 2));
		vect2.push_back(pow(f(i, G, r, sigma, T, x, K), 2));
	}	
    return mean(vect1) - mean(vect2);
}
	


double approx_1_derivee_variance(Vector Lambda, float r, float sigma, float T, float x, float K)
{
	// Lambda vecteur
	Vector G;
	double Y[2000];	
	double res = 0;
	
	for (int i = 0; i < 2000; i++)
	{
		G.push_back(gener_gauss());
	}
	
    for (int i = 0; i < 2000; i++)
	{
		Y[i] = (Lambda[i]-G[i])*exp(-Lambda[i]*G[i] + pow(Lambda[i],2)*0.5)*pow(f(i, G, r, sigma, T, x, K),2);
		res += Y[i];
	}
    return res/2000;
}


double approx_2_derivee_variance(Vector Lambda, float r, float sigma, float T, float x, float K)
{
	// Lambda vecteur
	Vector G;
	Vector G2;
	double Y[2000];
	double lam[2000]; 	
	double res = 0;
	
	for (int i = 0; i < 2000; i++)
	{	
		lam[i] = Lambda[i];
		G.push_back(gener_gauss());
		G2.push_back(G[i] + lam[i]);
	}
	
    for (int i = 0; i < 2000; i++)
	{
		Y[i] = -G[i]*exp(-2*Lambda[i]*G[i] - pow(Lambda[i],2))*pow(f(i, G2, r, sigma, T, x, K),2);
		res += Y[i];
	}
    return res/2000;
}


	
double F4(int n, double C, double beta, double g, double Lambda, float r, float sigma, float T, float x, float K)
{
	Vector G;
	for (int i = 0;i < n; i++) G.push_back(gener_gauss());
	
    return gamma(n,C,beta)*(Lambda-g)*exp(-Lambda*g+Lambda*Lambda*0.5)*pow(f(0, G, r, sigma, T, x, K), 2);
}

double F3(double Lambda, double g, float r, float sigma, float T, float x, float K)
{
	Vector G;
	for (int i = 0;i < 10; i++) G.push_back(gener_gauss());
	return (Lambda-g)*exp(-Lambda*g+Lambda*Lambda*0.5)*pow(f(0, G, r, sigma, T, x, K), 2);
}

void algorithme_rm_chen(double x_0, double MAX, int n, double C, double beta, float r, float sigma, float T, float x, float K,float (*F)(float, float))
{
	double X[n];
	double G[n];
	X[0] = x_0;
	
	for (int i = 0; i < n; i++)
	{
		G[i] = gener_gauss();
	}
	
	for (int i = 0; i < n-1; i++)
	{
		X[i+1] = X[i]-gamma(i,C,beta)*F(X[i],G[i+1], r, sigma, T, x, K);
		if ((max(-X[i+1], X[i+1])) >= MAX) X[i+1] = 0;
	}
	for (int i = 0; i < n; i++)
	{
		cout << X[i] << endl;
	}
}


	
	
	
int main()
{
	float r = 0.02;
	float sigma = 0.3;
	float T = 1;
	float x = 100;
	float K = 100;
	int x0 = 10;
	int n = 1000;
	float C = 1;
	float MAX = 5;
	float beta = 1;
	float G[n];
	for (int i = 0; i < n; i++)
	{
		G[i] = gener_gauss();
	}
	
	algorithme_rm(x0, n, C, beta,&F2);
	
	
}
	
	
	
	
	
	
	
	
