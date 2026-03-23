auto f = [=]( double t, double y, double& dy )
{
    dy = k * ( std::cos( t ) - y );
};

ponio::solve( f, ponio::runge_kutta::explicit_rkl1<5>(), y_0, t_span, dt, "ch_rkl1.txt"_fobs );
