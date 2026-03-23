auto f = [=]( double t, double y, double& dy )
{
    dy = k * ( std::cos( t ) - y );
};

ponio::solve( f, ponio::runge_kutta::rock::rock2(), y_0, t_span, dt, "ch_rock4.txt"_fobs );
