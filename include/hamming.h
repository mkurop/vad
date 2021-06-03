
    //! Generic fifth order symmetric cos window.
    /*!
     * \f$ w_i = a_0-a_1\ cos(2\pi i /(N-1))+a_2\ cos(4\pi i /(N-1))-a_3\ cos(6\pi i /(N-1))+a_4\ cos(8\pi i /(N-1))\f$
     *
     * \returns The cosinus window based on the <b>a</b> vector
     * \param[in] N Number of window taps
     * \param a A vector of cosinus coefficients
     */

    arma::Col<float> cos_win( const arma::uword N, const arma::Col<float> &a );

    //! Hamming window.
    /*! 
     * \f$ w_i = 0.54-0.46\ cos(2\pi i /(N-1))\f$
     * \param[in] N Nr of taps
     */
    arma::Col<float> hamming( const arma::uword N );
