auto f = [=]( double t, double y, double& dy )
{
    dy = k * ( std::cos( t ) - y );
};

ponio::solve( f, ponio::runge_kutta::rk_33(), y_0, t_span, dt, "ch_erk.txt"_fobs );
