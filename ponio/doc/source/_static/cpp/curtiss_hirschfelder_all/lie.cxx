auto f1 = [=]( double t, double y, double& dy )
{
    dy = k * std::cos( t );
};
auto f2 = [=]( double t, double y, double& dy )
{
    dy = -k * y;
};

auto pb = ponio::make_problem( f1, f2 );

auto lie = ponio::splitting::make_lie_tuple( std::make_pair( ponio::runge_kutta::rk_33(), 0.5 * dt ),
    std::make_pair( ponio::runge_kutta::rk_33(), 0.5 * dt ) );

ponio::solve( pb, lie, y_0, t_span, dt, "ch_lie.txt"_fobs );
