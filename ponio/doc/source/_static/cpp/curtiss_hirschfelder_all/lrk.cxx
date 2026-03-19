double L = -k;
auto N   = [=]( double t, double y, double& dy )
{
    dy = k * std::cos( t );
};

auto pb = ponio::make_lawson_problem( L, N );

auto my_exp = []( double x )
{
    return std::exp( x );
};

ponio::solve( pb, ponio::runge_kutta::lrk_33( my_exp ), y_0, t_span, dt, "ch_lrk.txt"_fobs );
