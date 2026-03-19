auto f_e = [=]( double t, double y, double& dy )
{
    dy = k * std::cos( t );
};
auto f_i = [=]( double t, double y, double& dy )
{
    dy = -k * y;
};
auto df_i = [=]( double t, double y )
{
    return -k;
};

auto pb = ponio::make_imex_jacobian_problem( f_e, f_i, df_i );

ponio::solve( pb, ponio::runge_kutta::imex_rk36_spi2(), y_0, t_span, dt, "ch_ark.txt"_fobs );
