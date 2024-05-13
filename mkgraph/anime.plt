# anime.plt
#-------------------------------------------------------------------------------
# プロット
#-------------------------------------------------------------------------------
do for[i=1000:100000:1000]{
    filename = sprintf("/data/sht/nakanog/taylor_test2/gnu_0/0_%d.d", i)
    #タイムステップを表示
    plot_title = sprintf('t*=%d',i)
    set title plot_title
    #plot[ｘ軸の範囲][y軸の範囲]
    splot [0:80][0:80] filename u 1:2:3
}