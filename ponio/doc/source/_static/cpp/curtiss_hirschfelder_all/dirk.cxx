auto f = [=]( double t, double y, double& dy )
{
    dy = k * ( std::cos( t ) - y );
};
auto df = [=]( double t, double /* y */ )
{
    return -k;
};

auto pb = ponio::make_implicit_problem( f, df );

ponio::solve( pb, ponio::runge_kutta::dirk34(), y_0, t_span, dt, "ch_dirk.txt"_fobs );
