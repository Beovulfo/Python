        function dTdZ(Z,gam,Zs,betha)
        real(8),intent(in):: Z,gam,Zs,betha
        real(8)::dTdZ
        dTdZ = gam/((1d0-Zs)*Zs)*0.5d0*(tanh(betha*(Zs-Z))+1d0)
        end function
        function Temper(H,Z,gam,Zs,betha,C)
        real(8),intent(in):: H,Z,gam,Zs,betha,C
        real(8)::Temper
        Temper =C+1d0+H+ (gam/((1d0-Zs)*Zs))/(2d0*betha)*(
     .  betha*Z-log(cosh(betha*(Zs-Z)))+log(cosh(betha*Zs)))
        end function
        function dZbdZ(Z,Zs,Zsbar,betha)
        real(8),intent(in):: Z,Zs,Zsbar,betha
        real(8)::dZbdZ
        dZbdZ = Zsbar/Zs+((1d0-Zsbar)/(1d0-Zs)-Zsbar/Zs)*
     .             0.5d0*(tanh(betha*(Z-Zs))+1d0)
        end function
        function dZbdZ2(Z,Zs,Zsbar,betha)
        real(8),intent(in):: Z,Zs,Zsbar,betha
        real(8)::dZbdZ2
        dZbdZ2 =((1d0-Zsbar)/(1d0-Zs)-Zsbar/Zs)*
     .              betha*0.5d0/(cosh(betha*(Z-Zs)))**2
        end function
