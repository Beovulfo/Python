        function dTdZ(Z,gam,Zs,Zsbar,betha)
        real(8),intent(in):: Z,gam,Zs,Zsbar,betha
        real(8)::dTdZ
        dTdZ = gam/((1d0-Zs)*Zsbar)*0.5d0*(tanh(betha*(Zsbar-Z))+1d0)
        !dTdZ = gam/((1d0-Zs)*Zs)*0.5d0*(dtanh(betha*(Zs-Z))+1d0)
        end function



        function Temper(H,Z,gam,Zs,Zsbar,betha,C)
        real(8),intent(in):: H,Z,gam,Zs,Zsbar,betha,C
        real(8)::Temper
        Temper = C+1d0 + H +0.5d0*(gam/((1-Zs)*Zsbar))*(Z +1/betha*(-log
     .   (cosh( betha*(Z-Zsbar)))+log(dcosh(betha*Zsbar))))
        end function
