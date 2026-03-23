auto f = [=]( double t, double y, double& dy )
{
    dy = k * ( std::cos( t ) - y );
};

ponio::solve( f, ponio::runge_kutta::rk54_6m().abs_tol( 1e-6 ).rel_tol( 1e-4 ), y_0, t_span, dt, "ch_dp.txt"_fobs );
