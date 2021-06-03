#include <hamming.h>

arma::Col<float> cos_win( const arma::uword N, const arma::Col<float> &a )
{
    arma::Col<float> h(N);
    for(arma::uword i=0; i<N; i++)
    {
        h[i] = a[0] - a[1]*std::cos(1.0*PI_2*i/(N-1)) + a[2]*std::cos(2.0*PI_2*i/(N-1)) \
               - a[3]*std::cos(3.0*PI_2*i/(N-1)) + a[4]*std::cos(4.0*PI_2*i/(N-1));
    }
    return h;
}

arma::Col<float> hamming( const arma::uword N )
{
    arma::vec a=arma::zeros<arma::Col<float>>(5);
    a[0] = 0.54;
    a[1] = 0.46;
    return cos_win(N,a);
}
