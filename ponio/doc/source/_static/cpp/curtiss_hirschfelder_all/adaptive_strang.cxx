auto f1 = [=]( double t, double y, double& dy )
{
    dy = k * std::cos( t );
};
auto f2 = [=]( double t, double y, double& dy )
{
    dy = -k * y;
};

auto pb = ponio::make_problem( f1, f2 );

auto strang = ponio::splitting::make_adaptive_strang_tuple( std::make_pair( ponio::runge_kutta::rk_33(), 0.5 * dt ),
    std::make_pair( ponio::runge_kutta::rk_33(), 0.5 * dt ) );

ponio::solve( pb, strang, y_0, t_span, dt, "ch_adaptive_strang.txt"_fobs );
