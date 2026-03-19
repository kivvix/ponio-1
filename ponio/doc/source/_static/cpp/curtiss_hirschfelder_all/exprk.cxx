double L = -k;
auto N   = [=]( double t, double y, double& dy )
{
    dy = k * std::cos( t );
};

auto pb = ponio::make_lawson_problem( L, N );

ponio::solve( pb, ponio::runge_kutta::exprk22(), y_0, t_span, dt, "ch_exprk.txt"_fobs );
