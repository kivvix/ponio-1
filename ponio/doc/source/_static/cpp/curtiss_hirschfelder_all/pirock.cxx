auto f1 = [=]( double t, double y, double& dy )
{
    dy = k * std::cos( t );
};
auto f2 = [=]( double t, double y, double& dy )
{
    dy = -k * y;
};
auto df2 = [=]( double t, double y )
{
    return -k;
};

auto pb = ponio::make_imex_jacobian_problem( f1, f2, df2 );

ponio::solve( pb, ponio::runge_kutta::pirock::pirock_a1(), y_0, t_span, dt, "ch_pirock.txt"_fobs );
