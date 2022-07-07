# define fixed axis-ranges
set xrange [-20:20]
set yrange [-10:10]
set zrange [-10:10]

set key noautotitle

# filename and n=number of lines of your data 
filedata = '../Outputs/PTMF_3D_10GeV.dat'
n = system(sprintf('cat %s | wc -l', filedata))

set term gif animate
set output '../Results/Trajectories10-99GeV.gif'

splot "../Outputs/BOX.dat" u 1:2:3 w l
do for [j=1:n] {
    set title 'time '.j
    replot'../Outputs/PTMF_3D_10GeV.dat' u 1:2:3 every ::j-1::j w p pt 7 ps 2, \
          '../Outputs/PTMF_3D_20GeV.dat' u 1:2:3 every ::1::j w l lw 2, \
          '../Outputs/PTMF_3D_20GeV.dat' u 1:2:3 every ::j-1::j w p pt 7 ps 2,\
          '../Outputs/PTMF_3D_30GeV.dat' u 1:2:3 every ::1::j w l lw 2, \
          '../Outputs/PTMF_3D_30GeV.dat' u 1:2:3 every ::j-1::j w p pt 7 ps 2, \
          '../Outputs/PTMF_3D_40GeV.dat' u 1:2:3 every ::1::j w l lw 2, \
          '../Outputs/PTMF_3D_40GeV.dat' u 1:2:3 every ::j-1::j w p pt 7 ps 2, \
          '../Outputs/PTMF_3D_50GeV.dat' u 1:2:3 every ::1::j w l lw 2, \
          '../Outputs/PTMF_3D_50GeV.dat' u 1:2:3 every ::j-1::j w p pt 7 ps 2, \
          '../Outputs/PTMF_3D_60GeV.dat' u 1:2:3 every ::1::j w l lw 2, \
          '../Outputs/PTMF_3D_60GeV.dat' u 1:2:3 every ::j-1::j w p pt 7 ps 2, \
          '../Outputs/PTMF_3D_70GeV.dat' u 1:2:3 every ::1::j w l lw 2, \
          '../Outputs/PTMF_3D_70GeV.dat' u 1:2:3 every ::j-1::j w p pt 7 ps 2, \
          '../Outputs/PTMF_3D_80GeV.dat' u 1:2:3 every ::1::j w l lw 2, \
          '../Outputs/PTMF_3D_80GeV.dat' u 1:2:3 every ::j-1::j w p pt 7 ps 2, \
          '../Outputs/PTMF_3D_90GeV.dat' u 1:2:3 every ::1::j w l lw 2, \
          '../Outputs/PTMF_3D_90GeV.dat' u 1:2:3 every ::j-1::j w p pt 7 ps 2, \
          '../Outputs/PTMF_3D_99GeV.dat' u 1:2:3 every ::1::j w l lw 2, \
          '../Outputs/PTMF_3D_99GeV.dat' u 1:2:3 every ::j-1::j w p pt 7 ps 2
    pause 0.02
}

# filename and n=number of lines of your data 
filedata = '../Outputs/PTMF_3D_10GeV.dat'
n = system(sprintf('cat %s | wc -l', filedata))

set term gif animate
set output '../Results/Trajectory10GeV.gif'

splot "../Outputs/BOX.dat" u 1:2:3 w l
do for [j=1:n] {
    set title 'time '.j
    replot '../Outputs/PTMF_3D_10GeV.dat' u 1:2:3 every ::1::j w l lw 2
    pause 0.02
}

# filename and n=number of lines of your data 
filedata = '../Outputs/PTMF_3D_30GeV.dat'
n = system(sprintf('cat %s | wc -l', filedata))

set term gif animate
set output '../Results/Trajectory20GeV.gif'

splot "../Outputs/BOX.dat" u 1:2:3 w l
do for [j=1:n] {
    set title 'time '.j
    replot '../Outputs/PTMF_3D_20GeV.dat' u 1:2:3 every ::1::j w l lw 2
    pause 0.02
}

# filename and n=number of lines of your data 
filedata = '../Outputs/PTMF_3D_50GeV.dat'
n = system(sprintf('cat %s | wc -l', filedata))

set term gif animate
set output '../Results/Trajectory50GeV.gif'

splot "../Outputs/BOX.dat" u 1:2:3 w l
do for [j=1:n] {
    set title 'time '.j
    replot '../Outputs/PTMF_3D_50GeV.dat' u 1:2:3 every ::1::j w l lw 2
    pause 0.02
}

# filename and n=number of lines of your data 
filedata = '../Outputs/PTMF_3D_70GeV.dat'
n = system(sprintf('cat %s | wc -l', filedata))

set term gif animate
set output '../Results/Trajectory70GeV.gif'

splot "../Outputs/BOX.dat" u 1:2:3 w l
do for [j=1:n] {
    set title 'time '.j
    replot '../Outputs/PTMF_3D_70GeV.dat' u 1:2:3 every ::1::j w l lw 2
    pause 0.02
}

# filename and n=number of lines of your data 
filedata = '../Outputs/PTMF_3D_90GeV.dat'
n = system(sprintf('cat %s | wc -l', filedata))

set term gif animate
set output '../Results/Trajectory90GeV.gif'

splot "../Outputs/BOX.dat" u 1:2:3 w l
do for [j=1:n] {
    set title 'time '.j
    replot '../Outputs/PTMF_3D_90GeV.dat' u 1:2:3 every ::1::j w l lw 2
    pause 0.02
}

# filename and n=number of lines of your data 
filedata = '../Outputs/PTMF_3D_99GeV.dat'
n = system(sprintf('cat %s | wc -l', filedata))

set term gif animate
set output '../Results/Trajectory99GeV.gif'

splot "../Outputs/BOX.dat" u 1:2:3 w l
do for [j=1:n] {
    set title 'time '.j
    replot '../Outputs/PTMF_3D_99GeV.dat' u 1:2:3 every ::1::j w l lw 2
    pause 0.02
}
