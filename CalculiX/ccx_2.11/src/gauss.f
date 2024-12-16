!
!     contains Gauss point information
!
!     gauss1d2: lin, 2-point integration (2 integration points)
!     gauss1d3: lin, 3-point integration (3 integration points)
!     gauss2d1: quad, 1-point integration (1 integration point)
!     gauss2d2: quad, 2-point integration (4 integration points)
!     gauss2d3: quad, 3-point integration (9 integration points)
!     gauss2d4: tri, 1 integration point
!     gauss2d5: tri, 3 integration points
!     gauss2d6: tri, 7 integration points
!     gauss3d1: hex, 1-point integration (1 integration point)
!     gauss3d2: hex, 2-point integration (8 integration points)
!     gauss3d3: hex, 3-point integration (27 integration points)
!     gauss3d4: tet, 1 integration point
!     gauss3d5: tet, 4 integration points
!     gauss3d6: tet, 15 integration points
!     gauss3d7: wedge, 2 integration points
!     gauss3d8: wedge, 9 integration points
!     gauss3d9: wedge, 18 integration points
!     gauss3d10: wedge, 6 integration points
!     gauss3d11: wedge, 1 integration points
!     gauss3d12: hex, 14 integration points (for c3d27)
!     gauss3d13: hex, 2x5x5=50 integration points (for beams)
!     gauss3d14: wedge, 1 integration point
!
!     weight2d1,... contains the weights
!
!
      real*8 gauss1d2(1,2),gauss1d3(1,3),
     &  gauss2d1(2,1),gauss2d2(2,4),gauss2d3(2,9),gauss2d4(2,1),
     &  gauss2d5(2,3),gauss3d1(3,1),gauss3d2(3,8),gauss3d3(3,27),
     &  gauss3d4(3,1),gauss3d5(3,4),gauss3d6(3,15),gauss3d7(3,2),
     &  gauss3d8(3,9),gauss3d9(3,18),gauss3d10(3,6),gauss3d11(3,1),
     &  gauss2d6(2,7),gauss3d12(3,14),gauss3d13(3,50),gauss3d14(3,1),
     &  weight1d2(2),weight1d3(3),weight2d1(1),weight2d2(4),
     &  weight2d3(9),weight2d4(1),weight2d5(3),weight3d1(1),
     &  weight3d2(8),weight3d3(27),weight3d4(1),weight3d5(4),
     &  weight3d6(15),weight3d7(2),weight3d8(9),weight3d9(18),
     &  weight3d10(6),weight3d11(1),weight2d6(7),weight3d12(14),
     &  weight3d13(50),weight3d14(1)
!
      gauss1d2=reshape(( /
     &  -0.577350269189626d0,0.577350269189626d0/),(/1,2/))
!
      gauss1d3=reshape(( /
     & -0.774596669241483d0,0.d0,0.774596669241483d0/),(/1,3/))
!
      gauss2d1=reshape((/0.,0./),(/2,1/))
!
!     the order of the Gauss points in gauss2d2 is important
!     and should not be changed (used to accelerate the code
!     for CAX8R axisymmetric elements in e_c3d_th.f)
!
      gauss2d2=reshape((/
     &  -0.577350269189626d0,-0.577350269189626d0,
     &  0.577350269189626d0,-0.577350269189626d0,
     &  -0.577350269189626d0,0.577350269189626d0,
     &  0.577350269189626d0,0.577350269189626d0/),(/2,4/))
!
      gauss2d3=reshape((/
     & -0.774596669241483d0,-0.774596669241483d0,
     & -0.d0,-0.774596669241483d0,
     & 0.774596669241483d0,-0.774596669241483d0,
     & -0.774596669241483d0,0.d0,
     & -0.d0,0.d0,
     & 0.774596669241483d0,0.d0,
     & -0.774596669241483d0,0.774596669241483d0,
     & -0.d0,0.774596669241483d0,
     & 0.774596669241483d0,0.774596669241483d0/),(/2,9/))
!
      gauss2d4=reshape((/0.333333333333333d0,0.333333333333333d0/),
     &  (/2,1/))
!
      gauss2d5=reshape((/
     & 0.166666666666667d0,0.166666666666667d0,
     & 0.666666666666667d0,0.166666666666667d0,
     & 0.166666666666667d0,0.666666666666667d0/),(/2,3/))
!
      gauss2d6=reshape((/
     & 0.333333333333333d0,0.333333333333333d0,
     & 0.797426985353087d0,0.101286507323456d0,
     & 0.101286507323456d0,0.797426985353087d0,
     & 0.101286507323456d0,0.101286507323456d0,
     & 0.470142064105115d0,0.059715871789770d0,
     & 0.059715871789770d0,0.470142064105115d0,
     & 0.470142064105115d0,0.470142064105115d0/),(/2,7/))
!
!
      gauss3d1=reshape((/0.,0.,0./),(/3,1/))
!
!     the order of the Gauss points in gauss3d2 is important
!     and should not be changed (used to accelerate the code
!     for CAX8R axisymmetric elements in e_c3d_th.f)
!
      gauss3d2=reshape((/
     &  -0.577350269189626d0,-0.577350269189626d0,-0.577350269189626d0,
     &  0.577350269189626d0,-0.577350269189626d0,-0.577350269189626d0,
     &  -0.577350269189626d0,0.577350269189626d0,-0.577350269189626d0,
     &  0.577350269189626d0,0.577350269189626d0,-0.577350269189626d0,
     &  -0.577350269189626d0,-0.577350269189626d0,0.577350269189626d0,
     &  0.577350269189626d0,-0.577350269189626d0,0.577350269189626d0,
     &  -0.577350269189626d0,0.577350269189626d0,0.577350269189626d0,
     &  0.577350269189626d0,0.577350269189626d0,0.577350269189626d0/),
     &  (/3,8/))
!
      gauss3d3=reshape((/
     & -0.774596669241483d0,-0.774596669241483d0,-0.774596669241483d0,
     & 0.d0,-0.774596669241483d0,-0.774596669241483d0,
     & 0.774596669241483d0,-0.774596669241483d0,-0.774596669241483d0,
     & -0.774596669241483d0,0.d0,-0.774596669241483d0,
     & 0.d0,0.d0,-0.774596669241483d0,
     & 0.774596669241483d0,0.d0,-0.774596669241483d0,
     & -0.774596669241483d0,0.774596669241483d0,-0.774596669241483d0,
     & 0.d0,0.774596669241483d0,-0.774596669241483d0,
     & 0.774596669241483d0,0.774596669241483d0,-0.774596669241483d0,
     & -0.774596669241483d0,-0.774596669241483d0,0.d0,
     & 0.d0,-0.774596669241483d0,0.d0,
     & 0.774596669241483d0,-0.774596669241483d0,0.d0,
     & -0.774596669241483d0,0.d0,0.d0,
     & 0.d0,0.d0,0.d0,
     & 0.774596669241483d0,0.d0,0.d0,
     & -0.774596669241483d0,0.774596669241483d0,0.d0,
     & 0.d0,0.774596669241483d0,0.d0,
     & 0.774596669241483d0,0.774596669241483d0,0.d0,
     & -0.774596669241483d0,-0.774596669241483d0,0.774596669241483d0,
     & 0.d0,-0.774596669241483d0,0.774596669241483d0,
     & 0.774596669241483d0,-0.774596669241483d0,0.774596669241483d0,
     & -0.774596669241483d0,0.d0,0.774596669241483d0,
     & 0.d0,0.d0,0.774596669241483d0,
     & 0.774596669241483d0,0.d0,0.774596669241483d0,
     & -0.774596669241483d0,0.774596669241483d0,0.774596669241483d0,
     & 0.d0,0.774596669241483d0,0.774596669241483d0,
     & 0.774596669241483d0,0.774596669241483d0,0.774596669241483d0/),
     &  (/3,27/))
!
      gauss3d4=reshape((/0.25d0,0.25d0,0.25d0/),(/3,1/))
!
      gauss3d5=reshape((/
     & 0.138196601125011d0,0.138196601125011d0,0.138196601125011d0,
     & 0.585410196624968d0,0.138196601125011d0,0.138196601125011d0,
     & 0.138196601125011d0,0.585410196624968d0,0.138196601125011d0,
     & 0.138196601125011d0,0.138196601125011d0,0.585410196624968d0/),
     &  (/3,4/))
!
      gauss3d6=reshape((/
     & 0.25d0,0.25d0,0.25d0,
     & 0.091971078052723d0,0.091971078052723d0,0.091971078052723d0,
     & 0.724086765841831d0,0.091971078052723d0,0.091971078052723d0,
     & 0.091971078052723d0,0.724086765841831d0,0.091971078052723d0,
     & 0.091971078052723d0,0.091971078052723d0,0.724086765841831d0,
     & 0.319793627829630d0,0.319793627829630d0,0.319793627829630d0,
     & 0.040619116511110d0,0.319793627829630d0,0.319793627829630d0,
     & 0.319793627829630d0,0.040619116511110d0,0.319793627829630d0,
     & 0.319793627829630d0,0.319793627829630d0,0.040619116511110d0,
     & 0.056350832689629d0,0.056350832689629d0,0.443649167310371d0,
     & 0.443649167310371d0,0.056350832689629d0,0.056350832689629d0,
     & 0.443649167310371d0,0.443649167310371d0,0.056350832689629d0,
     & 0.056350832689629d0,0.443649167310371d0,0.443649167310371d0,
     & 0.056350832689629d0,0.443649167310371d0,0.056350832689629d0,
     & 0.443649167310371d0,0.056350832689629d0,0.443649167310371d0/),
     &  (/3,15/))
!
      gauss3d7=reshape((/
     & 0.333333333333333d0,0.333333333333333d0,-0.577350269189626d0,
     & 0.333333333333333d0,0.333333333333333d0,0.577350269189626d0/),
     &  (/3,2/))
!
      gauss3d8=reshape((/
     & 0.166666666666667d0,0.166666666666667d0,-0.774596669241483d0,
     & 0.666666666666667d0,0.166666666666667d0,-0.774596669241483d0,
     & 0.166666666666667d0,0.666666666666667d0,-0.774596669241483d0,
     & 0.166666666666667d0,0.166666666666667d0,0.d0,
     & 0.666666666666667d0,0.166666666666667d0,0.d0,
     & 0.166666666666667d0,0.666666666666667d0,0.d0,
     & 0.166666666666667d0,0.166666666666667d0,0.774596669241483d0,
     & 0.666666666666667d0,0.166666666666667d0,0.774596669241483d0,
     & 0.166666666666667d0,0.666666666666667d0,0.774596669241483d0/),
     &  (/3,9/))
!
      gauss3d9=reshape((/
     & 0.166666666666667d0,0.166666666666667d0,-0.774596669241483d0,
     & 0.166666666666667d0,0.666666666666667d0,-0.774596669241483d0,
     & 0.666666666666667d0,0.166666666666667d0,-0.774596669241483d0,
     & 0.000000000000000d0,0.500000000000000d0,-0.774596669241483d0,
     & 0.500000000000000d0,0.000000000000000d0,-0.774596669241483d0,
     & 0.500000000000000d0,0.500000000000000d0,-0.774596669241483d0,
     & 0.166666666666667d0,0.166666666666667d0,0.d0,
     & 0.166666666666667d0,0.666666666666667d0,0.d0,
     & 0.666666666666667d0,0.166666666666667d0,0.d0,
     & 0.000000000000000d0,0.500000000000000d0,0.d0,
     & 0.500000000000000d0,0.000000000000000d0,0.d0,
     & 0.500000000000000d0,0.500000000000000d0,0.d0,
     & 0.166666666666667d0,0.166666666666667d0,0.774596669241483d0,
     & 0.166666666666667d0,0.666666666666667d0,0.774596669241483d0,
     & 0.666666666666667d0,0.166666666666667d0,0.774596669241483d0,
     & 0.000000000000000d0,0.500000000000000d0,0.774596669241483d0,
     & 0.500000000000000d0,0.000000000000000d0,0.774596669241483d0,
     & 0.500000000000000d0,0.500000000000000d0,0.774596669241483d0/),
     &  (/3,18/))
!
      gauss3d10=reshape((/
     & 0.166666666666667d0,0.166666666666667d0,-0.577350269189626d0,
     & 0.666666666666667d0,0.166666666666667d0,-0.577350269189626d0,
     & 0.166666666666667d0,0.666666666666667d0,-0.577350269189626d0,
     & 0.166666666666667d0,0.166666666666667d0,0.577350269189626d0,
     & 0.666666666666667d0,0.166666666666667d0,0.577350269189626d0,
     & 0.166666666666667d0,0.666666666666667d0,0.577350269189626d0/),
     &  (/3,6/))
!
      gauss3d11=reshape((/
     & 0.333333333333333d0,0.333333333333333d0,0.d0/),(/3,1/))
!
      gauss3d12=reshape((/
     & -0.758786910639328d0,-0.758786910639328d0,-0.758786910639328d0,
     & 0.758786910639328d0,-0.758786910639328d0,-0.758786910639328d0,
     & 0.758786910639328d0,0.758786910639328d0,-0.758786910639328d0,
     & -0.758786910639328d0,0.758786910639328d0,-0.758786910639328d0,
     & -0.758786910639328d0,-0.758786910639328d0,0.758786910639328d0,
     & 0.758786910639328d0,-0.758786910639328d0,0.758786910639328d0,
     & 0.758786910639328d0,0.758786910639328d0,0.758786910639328d0,
     & -0.758786910639328d0,0.758786910639328d0,0.758786910639328d0,
     & 0.d0,-0.795822425754222d0,0.d0,
     & 0.795822425754222d0,0.d0,0.d0,
     & 0.d0,0.795822425754222d0,0.d0,
     & -0.795822425754222d0,0.d0,0.d0,
     & 0.d0,0.d0,-0.795822425754222d0,
     & 0.d0,0.d0,0.795822425754222d0/),
     &  (/3,14/))
!
      gauss3d13=reshape(( /
     &-0.577350269189626d+00,-0.5773502692d+00,-0.5773502692d+00,
     & 0.577350269189626d+00,-0.5773502692d+00,-0.5773502692d+00,
     &-0.577350269189626d+00, 0.5773502692d+00,-0.5773502692d+00,
     & 0.577350269189626d+00, 0.5773502692d+00,-0.5773502692d+00,
     &-0.577350269189626d+00,-0.5773502692d+00, 0.5773502692d+00,
     & 0.577350269189626d+00,-0.5773502692d+00, 0.5773502692d+00,
     &-0.577350269189626d+00, 0.5773502692d+00, 0.5773502692d+00,
     & 0.577350269189626d+00, 0.5773502692d+00, 0.5773502692d+00,
     &-0.577350269189626d+00,-0.9258200998d+00,-0.9258200998d+00,
     &-0.577350269189626d+00,-0.5773502692d+00,-0.9258200998d+00,
     &-0.577350269189626d+00, 0.0000000000d+00,-0.9258200998d+00,
     &-0.577350269189626d+00, 0.5773502692d+00,-0.9258200998d+00,
     &-0.577350269189626d+00, 0.9258200998d+00,-0.9258200998d+00,
     &-0.577350269189626d+00,-0.9258200998d+00,-0.5773502692d+00,
     &-0.577350269189626d+00, 0.0000000000d+00,-0.5773502692d+00,
     &-0.577350269189626d+00, 0.9258200998d+00,-0.5773502692d+00,
     &-0.577350269189626d+00,-0.9258200998d+00, 0.0000000000d+00,
     &-0.577350269189626d+00,-0.5773502692d+00, 0.0000000000d+00,
     &-0.577350269189626d+00, 0.0000000000d+00, 0.0000000000d+00,
     &-0.577350269189626d+00, 0.5773502692d+00, 0.0000000000d+00,
     &-0.577350269189626d+00, 0.9258200998d+00, 0.0000000000d+00,
     &-0.577350269189626d+00,-0.9258200998d+00, 0.5773502692d+00,
     &-0.577350269189626d+00, 0.0000000000d+00, 0.5773502692d+00,
     &-0.577350269189626d+00, 0.9258200998d+00, 0.5773502692d+00,
     &-0.577350269189626d+00,-0.9258200998d+00, 0.9258200998d+00,
     &-0.577350269189626d+00,-0.5773502692d+00, 0.9258200998d+00,
     &-0.577350269189626d+00, 0.0000000000d+00, 0.9258200998d+00,
     &-0.577350269189626d+00, 0.5773502692d+00, 0.9258200998d+00,
     &-0.577350269189626d+00, 0.9258200998d+00, 0.9258200998d+00,
     & 0.577350269189626d+00,-0.9258200998d+00,-0.9258200998d+00,
     & 0.577350269189626d+00,-0.5773502692d+00,-0.9258200998d+00,
     & 0.577350269189626d+00, 0.0000000000d+00,-0.9258200998d+00,
     & 0.577350269189626d+00, 0.5773502692d+00,-0.9258200998d+00,
     & 0.577350269189626d+00, 0.9258200998d+00,-0.9258200998d+00,
     & 0.577350269189626d+00,-0.9258200998d+00,-0.5773502692d+00,
     & 0.577350269189626d+00, 0.0000000000d+00,-0.5773502692d+00,
     & 0.577350269189626d+00, 0.9258200998d+00,-0.5773502692d+00,
     & 0.577350269189626d+00,-0.9258200998d+00, 0.0000000000d+00,
     & 0.577350269189626d+00,-0.5773502692d+00, 0.0000000000d+00,
     & 0.577350269189626d+00, 0.0000000000d+00, 0.0000000000d+00,
     & 0.577350269189626d+00, 0.5773502692d+00, 0.0000000000d+00,
     & 0.577350269189626d+00, 0.9258200998d+00, 0.0000000000d+00,
     & 0.577350269189626d+00,-0.9258200998d+00, 0.5773502692d+00,
     & 0.577350269189626d+00, 0.0000000000d+00, 0.5773502692d+00,
     & 0.577350269189626d+00, 0.9258200998d+00, 0.5773502692d+00,
     & 0.577350269189626d+00,-0.9258200998d+00, 0.9258200998d+00,
     & 0.577350269189626d+00,-0.5773502692d+00, 0.9258200998d+00,
     & 0.577350269189626d+00, 0.0000000000d+00, 0.9258200998d+00,
     & 0.577350269189626d+00, 0.5773502692d+00, 0.9258200998d+00,
     & 0.577350269189626d+00, 0.9258200998d+00, 0.9258200998d+00/),
     & (/3,50/))
!
      gauss3d14=reshape((/
     & 0.333333333333333d0,0.333333333333333d0,0.d0/),
     &  (/3,1/))
!
      weight1d2=(/1.d0,1.d0/)
!
      weight1d3=(/0.555555555555555d0,0.888888888888888d0,
     &  0.555555555555555d0/)
!
      weight2d1=(/4.d0/)
!
      weight2d2=(/1.d0,1.d0,1.d0,1.d0/)
!
      weight2d3=(/
     &  0.308641975308642d0,0.493827160493827d0,0.308641975308642d0,
     &  0.493827160493827d0,0.790123456790123d0,0.493827160493827d0,
     &  0.308641975308642d0,0.493827160493827d0,0.308641975308642d0/)
!
      weight2d4=(/0.5d0/)
!
      weight2d5=(/
     &  0.166666666666666d0,0.166666666666666d0,0.166666666666666d0/)
!
c      weight2d6=(/
c     & 0.225000000000000d0,0.125939180544827d0,0.125939180544827d0,
c     & 0.125939180544827d0,0.132394152788506d0,0.132394152788506d0,
c     & 0.132394152788506d0/)
!
      weight2d6=(/
     & 0.112500000000000d0,0.062969590272413d0,0.062969590272413d0,
     & 0.062969590272413d0,0.066197076394253d0,0.066197076394253d0,
     & 0.066197076394253d0/)
!
      weight3d1=(/8.d0/)
!
      weight3d2=(/1.d0,1.d0,1.d0,1.d0,1.d0,1.d0,1.d0,1.d0/)
!
      weight3d3=(/
     &  0.171467764060357d0,0.274348422496571d0,0.171467764060357d0,
     &  0.274348422496571d0,0.438957475994513d0,0.274348422496571d0,
     &  0.171467764060357d0,0.274348422496571d0,0.171467764060357d0,
     &  0.274348422496571d0,0.438957475994513d0,0.274348422496571d0,
     &  0.438957475994513d0,0.702331961591221d0,0.438957475994513d0,
     &  0.274348422496571d0,0.438957475994513d0,0.274348422496571d0,
     &  0.171467764060357d0,0.274348422496571d0,0.171467764060357d0,
     &  0.274348422496571d0,0.438957475994513d0,0.274348422496571d0,
     &  0.171467764060357d0,0.274348422496571d0,0.171467764060357d0/)
!
      weight3d4=(/0.166666666666667d0/)
!
      weight3d5=(/
     &  0.041666666666667d0,0.041666666666667d0,0.041666666666667d0,
     &  0.041666666666667d0/)
!
      weight3d6=(/
     &  0.019753086419753d0,0.011989513963170d0,0.011989513963170d0,
     &  0.011989513963170d0,0.011989513963170d0,0.011511367871045d0,
     &  0.011511367871045d0,0.011511367871045d0,0.011511367871045d0,
     &  0.008818342151675d0,0.008818342151675d0,0.008818342151675d0,
     &  0.008818342151675d0,0.008818342151675d0,0.008818342151675d0/)
!
      weight3d7=(/0.5d0,0.5d0/)
!
      weight3d8=(/
     &  0.092592592592593d0,0.092592592592593d0,0.092592592592593d0,
     &  0.148148148148148d0,0.148148148148148d0,0.148148148148148d0,
     &  0.092592592592593d0,0.092592592592593d0,0.092592592592593d0/)
!
      weight3d9=(/
     &  0.083333333333333d0,0.083333333333333d0,0.083333333333333d0,
     &  0.009259259259259d0,0.009259259259259d0,0.009259259259259d0,
     &  0.133333333333333d0,0.133333333333333d0,0.133333333333333d0,
     &  0.014814814814815d0,0.014814814814815d0,0.014814814814815d0,
     &  0.083333333333333d0,0.083333333333333d0,0.083333333333333d0,
     &  0.009259259259259d0,0.009259259259259d0,0.009259259259259d0/)
!
      weight3d10=(/
     &  0.166666666666666d0,0.166666666666666d0,0.166666666666666d0,
     &  0.166666666666666d0,0.166666666666666d0,0.166666666666666d0/)
!
      weight3d11=(/1.d0/)
!
      weight3d12=(/
     &  0.335180055401662d0,0.335180055401662d0,0.335180055401662d0,
     &  0.335180055401662d0,0.335180055401662d0,0.335180055401662d0,
     &  0.335180055401662d0,0.335180055401662d0,0.886426592797794d0,
     &  0.886426592797794d0,0.886426592797794d0,0.886426592797794d0,
     &  0.886426592797794d0,0.886426592797794d0/)
!
      weight3d13=(/
     & 0.240991735537190d+00, 0.240991735537190d+00,
     & 0.240991735537190d+00, 0.240991735537190d+00,
     & 0.240991735537190d+00, 0.240991735537190d+00,
     & 0.240991735537190d+00, 0.240991735537190d+00,
     & 0.391960004081216d-01, 0.971900826446281d-01,
     & 0.123187429854097d+00, 0.971900826446281d-01,
     & 0.391960004081216d-01, 0.971900826446281d-01,
     & 0.305454545454545d+00,
     & 0.971900826446281d-01,
     & 0.123187429854097d+00, 0.305454545454545d+00,
     & 0.387160493827161d+00, 0.305454545454545d+00,
     & 0.123187429854097d+00, 0.971900826446281d-01,
     & 0.305454545454545d+00,
     & 0.971900826446281d-01,
     & 0.391960004081216d-01, 0.971900826446281d-01,
     & 0.123187429854097d+00, 0.971900826446281d-01,
     & 0.391960004081216d-01, 0.391960004081216d-01,
     & 0.971900826446281d-01, 0.123187429854097d+00,
     & 0.971900826446281d-01, 0.391960004081216d-01,
     & 0.971900826446281d-01, 
     & 0.305454545454545d+00, 
     & 0.971900826446281d-01, 0.123187429854097d+00,
     & 0.305454545454545d+00, 0.387160493827161d+00,
     & 0.305454545454545d+00, 0.123187429854097d+00,
     & 0.971900826446281d-01, 
     & 0.305454545454545d+00, 
     & 0.971900826446281d-01, 0.391960004081216d-01,
     & 0.971900826446281d-01, 0.123187429854097d+00,
     & 0.971900826446281d-01, 0.391960004081216d-01/)
!
      weight3d14=(/1.d0/)
!
