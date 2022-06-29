# define fixed axis-ranges
set xrange [-20:20]
set yrange [-10:10]
set zrange [-10:10]

# filename and n=number of lines of your data 
filedata = 'PTMF_3D_10GeV.dat'
n = system(sprintf('cat %s | wc -l', filedata))

do for [j=1:n] {
    set title 'time '.j
    splot 'PTMF_3D_10GeV.dat' u 1:2:3 every ::1::j w l lw 2, \
          'PTMF_3D_10GeV.dat' u 1:2:3 every ::j-1::j w p pt 7 ps 2, \
          'PTMF_3D_20GeV.dat' u 1:2:3 every ::1::j w l lw 2, \
          'PTMF_3D_20GeV.dat' u 1:2:3 every ::j-1::j w p pt 7 ps 2,\
          'PTMF_3D_30GeV.dat' u 1:2:3 every ::1::j w l lw 2, \
          'PTMF_3D_30GeV.dat' u 1:2:3 every ::j-1::j w p pt 7 ps 2, \
          'PTMF_3D_40GeV.dat' u 1:2:3 every ::1::j w l lw 2, \
          'PTMF_3D_40GeV.dat' u 1:2:3 every ::j-1::j w p pt 7 ps 2, \
          'PTMF_3D_50GeV.dat' u 1:2:3 every ::1::j w l lw 2, \
          'PTMF_3D_50GeV.dat' u 1:2:3 every ::j-1::j w p pt 7 ps 2, \
          'PTMF_3D_60GeV.dat' u 1:2:3 every ::1::j w l lw 2, \
          'PTMF_3D_60GeV.dat' u 1:2:3 every ::j-1::j w p pt 7 ps 2, \
          'PTMF_3D_70GeV.dat' u 1:2:3 every ::1::j w l lw 2, \
          'PTMF_3D_70GeV.dat' u 1:2:3 every ::j-1::j w p pt 7 ps 2, \
          'PTMF_3D_80GeV.dat' u 1:2:3 every ::1::j w l lw 2, \
          'PTMF_3D_80GeV.dat' u 1:2:3 every ::j-1::j w p pt 7 ps 2, \
          'PTMF_3D_90GeV.dat' u 1:2:3 every ::1::j w l lw 2, \
          'PTMF_3D_90GeV.dat' u 1:2:3 every ::j-1::j w p pt 7 ps 2, \
          'PTMF_3D_99GeV.dat' u 1:2:3 every ::1::j w l lw 2, \
          'PTMF_3D_99GeV.dat' u 1:2:3 every ::j-1::j w p pt 7 ps 2
    pause 0.02
}
