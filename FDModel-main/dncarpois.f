C*****************
      subroutine pois (iflg,nperod,n,mperod,m,a,b,c,idimy,y,w)
      implicit real*8 (A-H,O-Z)
      external        trid       ,tridp
      dimension       y(idimy,*)
      dimension       w(*)       ,b(*)       ,a(*)       ,c(*)
      iwdim1 = 5.25D0*n+1
      iwdim2 = iwdim1+m
      iwdim3 = iwdim2+m
      iwdim4 = iwdim3+m
      iwdim5 = iwdim4+m
      do 102 i=1,m
         a(i) = -a(i)
         c(i) = -c(i)
  101    w(iwdim5+i-1) = -b(i)+2.D0
  102 continue
      if (mperod .eq. 0) go to 103
      call poisgn (nperod,n,iflg,m,a,w(iwdim5),c,idimy,y,w(1),
     1             w(iwdim1),w(iwdim2),w(iwdim3),w(iwdim4),trid)
      go to 104
  103 call poisgn (nperod,n,iflg,m,a,w(iwdim5),c,idimy,y,w(1),
     1             w(iwdim1),w(iwdim2),w(iwdim3),w(iwdim4),tridp)
  104 do 106 i=1,m
         a(i) = -a(i)
  105    c(i) = -c(i)
  106 continue
  107 return
      end
C**************************
      subroutine poisgn (nperod,n,iflg,m,ba,bb,bc,idimq,q,twocos,b,t,d,
     1                   w,tri)
      dimension       q(idimq,*)
      dimension       ba(*)      ,bb(*)      ,bc(*)      ,twocos(*)  ,
     1                b(*)       ,t(*)       ,d(*)       ,w(*)
      if (iflg .ne. 0) go to 101
c
c     do initialization if first time through.
c
      call poinit (nperod,n,iex2,iex3,iex5,n2pw,n2p3p,n2p3p5,ks3,ks5,
     1             npstop,twocos)
  101 mm1 = m-1
      do 104 j=1,n
         do 102 i=1,m
            b(i) = -q(i,j)
  102    continue
         call tri (1,1,1,m,mm1,ba,bb,bc,b,twocos,d,w)
         do 103 i=1,m
            q(i,j) = b(i)
  103    continue
  104 continue
      np = 1
c
c     start reduction for powers of two
c
      if (iex2 .eq. 0) go to 112
      l = iex2
      if (iex3+iex5 .ne. 0) go to 105
      if (nperod .eq. 1) l = l-1
      if (l .eq. 0) go to 138
  105 do 111 kount=1,l
         k = np
         np = 2*np
         k2 = np-1
         k3 = np
         k4 = 2*np-1
         jstart = np
         if (nperod .eq. 3) jstart = 1
         do 110 j=jstart,n,np
            jm1 = j-k
            jp1 = j+k
            if (j .eq. 1) jm1 = jp1
            if (j .ne. n) go to 106
            jp1 = jm1
            if (nperod .eq. 0) jp1 = k
  106       do 107 i=1,m
               b(i) = q(i,jm1)+q(i,jp1)
  107       continue
            call tri (k,k2,1,m,mm1,ba,bb,bc,b,twocos,d,w)
            do 108 i=1,m
               t(i) = q(i,j)+b(i)
               b(i) = t(i)
  108       continue
            call tri (k3,k4,1,m,mm1,ba,bb,bc,b,twocos,d,w)
            do 109 i=1,m
               q(i,j) = t(i)+2.*b(i)
  109       continue
  110    continue
  111 continue
c
c     start reduction for powers of three
c
  112 if (iex3 .eq. 0) go to 124
      l = iex3
      if (iex5 .ne. 0) go to 113
      if (nperod .eq. 1) l = l-1
      if (l .eq. 0) go to 138
  113 k2 = np-1
      do 123 kount=1,l
         k = np
         np = 3*np
         k1 = k2+1
         k2 = k2+k
         k3 = k2+1
         k4 = k2+np
         jstart = np
         if (nperod .eq. 3) jstart = 1
         do 122 j=jstart,n,np
            if (j .ne. 1) go to 114
            jm1 = j+k
            jm2 = jm1+k
            go to 116
  114       jm1 = j-k
            jm2 = jm1-k
            if (j .ne. n) go to 116
            if (nperod .eq. 0) go to 115
            jp1 = jm1
            jp2 = jm2
            go to 117
  115       jp1 = k
            jp2 = jp1+k
            go to 117
  116       jp1 = j+k
            jp2 = jp1+k
  117       do 118 i=1,m
               b(i) = 2.*q(i,j)+q(i,jm2)+q(i,jp2)
  118       continue
            call tri (k1,k2,1,m,mm1,ba,bb,bc,b,twocos,d,w)
            do 119 i=1,m
               t(i) = b(i)+q(i,jm1)+q(i,jp1)
               b(i) = t(i)
  119       continue
            call tri (k1,k2,1,m,mm1,ba,bb,bc,b,twocos,d,w)
            do 120 i=1,m
               q(i,j) = q(i,j)+b(i)
               b(i) = t(i)
  120       continue
            call tri (k3,k4,1,m,mm1,ba,bb,bc,b,twocos,d,w)
            do 121 i=1,m
               q(i,j) = q(i,j)+3.*b(i)
  121       continue
  122    continue
  123 continue
c
c     start reduction for powers of five
c
  124 l = iex5
      if (nperod .eq. 1) l = l-1
      if (l .le. 0) go to 138
      k2 = (n2pw+n2p3p)/2-1
      do 137 kount=1,l
         k = np
         np = 5*np
         k1 = k2+1
         k2 = k2+k
         k3 = k2+1
         k4 = k2+np
         jstart = np
         if (nperod .eq. 3) jstart = 1
         do 136 j=jstart,n,np
            if (j .ne. 1) go to 125
            jm1 = j+k
            jm2 = jm1+k
            jm3 = jm2+k
            jm4 = jm3+k
            go to 127
  125       jm1 = j-k
            jm2 = jm1-k
            jm3 = jm2-k
            jm4 = jm3-k
            if (j .ne. n) go to 127
            if (nperod .eq. 0) go to 126
            jp1 = jm1
            jp2 = jm2
            jp3 = jm3
            jp4 = jm4
            go to 128
  126       jp1 = k
            jp2 = jp1+k
            jp3 = jp2+k
            jp4 = jp3+k
            go to 128
  127       jp1 = j+k
            jp2 = jp1+k
            jp3 = jp2+k
            jp4 = jp3+k
  128       do 129 i=1,m
               b(i) = 6.*q(i,j)+4.*(q(i,jm2)+q(i,jp2))+q(i,jm4)+q(i,jp4)
  129       continue
            call tri (k1,k2,1,m,mm1,ba,bb,bc,b,twocos,d,w)
            do 130 i=1,m
               b(i) = b(i)+3.*(q(i,jm1)+q(i,jp1))+q(i,jm3)+q(i,jp3)
  130       continue
            call tri (k1,k2,1,m,mm1,ba,bb,bc,b,twocos,d,w)
            do 131 i=1,m
               t(i) = b(i)
               b(i) = 2.*q(i,j)+q(i,jm2)+q(i,jp2)+b(i)
  131       continue
            call tri (k1,k2,1,m,mm1,ba,bb,bc,b,twocos,d,w)
            do 132 i=1,m
               b(i) = b(i)+q(i,jm1)+q(i,jp1)
  132       continue
            call tri (k1,k2,1,m,mm1,ba,bb,bc,b,twocos,d,w)
            do 133 i=1,m
               temp = b(i)+q(i,j)
               b(i) = 4.*q(i,j)+3.*(q(i,jm2)+q(i,jp2))+q(i,jm4)+
     1                q(i,jp4)-t(i)
               q(i,j) = temp
  133       continue
            call tri (k1,k2,1,m,mm1,ba,bb,bc,b,twocos,d,w)
            do 134 i=1,m
               b(i) = b(i)+2.*(q(i,jm1)+q(i,jp1))+q(i,jm3)+q(i,jp3)
  134       continue
            call tri (k3,k4,1,m,mm1,ba,bb,bc,b,twocos,d,w)
            do 135 i=1,m
               q(i,j) = q(i,j)+5.*b(i)
  135       continue
  136    continue
  137 continue
c
c     initial phase of back substitution
c
  138 if (nperod.eq.1 .or. nperod.eq.2) go to 147
      if (nperod .eq. 0) go to 141
      do 139 i=1,m
         b(i) = 2.*q(i,1)
  139 continue
      call tri (k3,k4,1,m,mm1,ba,bb,bc,b,twocos,d,w)
      do 140 i=1,m
         q(i,n) = q(i,n)+b(i)
         b(i) = 4.*q(i,n)
  140 continue
      call tri (k4+1,k4+2*np-1,2,m,mm1,ba,bb,bc,b,twocos,d,w)
      go to 143
  141 do 142 i=1,m
         b(i) = 2.*q(i,n)
  142 continue
      call tri (k4+1,k4+np-1,2,m,mm1,ba,bb,bc,b,twocos,d,w)
  143 do 144 i=1,m
         q(i,n) = q(i,n)+b(i)
  144 continue
      if (nperod .eq. 0) go to 147
      do 145 i=1,m
         b(i) = 2.*q(i,n)
  145 continue
      call tri (k3,k4,1,m,mm1,ba,bb,bc,b,twocos,d,w)
      do 146 i=1,m
         q(i,1) = q(i,1)+b(i)
  146 continue
c
c     regular back substitution for powers of five
c
  147 continue
      if (iex5 .eq. 0) go to 171
      np = n2p3p5
      k8 = ks5
      if (nperod .eq. 1) k3 = k3+np/5
      do 170 kount=1,iex5
         k = np
         np = np/5
         k4 = k3-1
         k3 = k3-np
         k1 = k8+1
         k2 = k8+2*np
         k5 = k2+1
         k6 = k2+4*np
         k7 = k6+1
         k8 = k6+2*np
         jstart = np
         if (nperod .eq. 3) jstart = 1+np
         do 169 j=jstart,n,k
            jm1 = j-np
            jp1 = j+np
            jp2 = jp1+np
            jp3 = jp2+np
            jp4 = jp3+np
            if (jm1 .ne. 0) go to 151
            if (nperod .eq. 0) go to 149
            do 148 i=1,m
               b(i) = q(i,jp1)
  148       continue
            go to 153
  149       do 150 i=1,m
               b(i) = q(i,jp1)+q(i,n)
  150       continue
            go to 153
  151       do 152 i=1,m
               b(i) = q(i,jp1)+q(i,jm1)
  152       continue
  153       call tri (k3,k4,1,m,mm1,ba,bb,bc,b,twocos,d,w)
            do 154 i=1,m
               q(i,j) = q(i,j)+b(i)
               b(i) = q(i,jp1)+q(i,jp3)
  154       continue
            call tri (k3,k4,1,m,mm1,ba,bb,bc,b,twocos,d,w)
            do 155 i=1,m
               q(i,jp2) = q(i,jp2)+b(i)
  155       continue
            if (jp4 .gt. n) go to 157
            do 156 i=1,m
               b(i) = 2.*q(i,jp2)+q(i,j)+q(i,jp4)
  156       continue
            go to 159
  157       do 158 i=1,m
               b(i) = 2.*q(i,jp2)+q(i,j)
  158       continue
  159       call tri (k1,k2,1,m,mm1,ba,bb,bc,b,twocos,d,w)
            do 160 i=1,m
               q(i,jp2) = q(i,jp2)+b(i)
               b(i) = q(i,jp2)+q(i,j)
  160       continue
            call tri (k5,k6,1,m,mm1,ba,bb,bc,b,twocos,d,w)
            do 161 i=1,m
               q(i,jp2) = q(i,jp2)+b(i)
               b(i) = q(i,j)+q(i,jp2)
  161       continue
            call tri (k7,k8,1,m,mm1,ba,bb,bc,b,twocos,d,w)
            do 162 i=1,m
               q(i,j) = q(i,j)+b(i)
               b(i) = q(i,j)+q(i,jp2)
  162       continue
            call tri (k3,k4,1,m,mm1,ba,bb,bc,b,twocos,d,w)
            do 163 i=1,m
               q(i,jp1) = q(i,jp1)+b(i)
  163       continue
            if (jp4 .gt. n) go to 165
            do 164 i=1,m
               b(i) = q(i,jp2)+q(i,jp4)
  164       continue
            go to 167
  165       do 166 i=1,m
               b(i) = q(i,jp2)
  166       continue
  167       call tri (k3,k4,1,m,mm1,ba,bb,bc,b,twocos,d,w)
            do 168 i=1,m
               q(i,jp3) = q(i,jp3)+b(i)
  168       continue
  169    continue
  170 continue
c
c     regular back substitution for powers of three
c
  171 if (iex3 .eq. 0) go to 191
      np = n2p3p
      k2 = ks3
      if (nperod.eq.1 .and. iex5.eq.0) k3 = k3+np/3
      do 190 kount=1,iex3
         k = np
         np = np/3
         k4 = k3-1
         k3 = k3-np
         k1 = k2+1
         k2 = k2+2*np
         jstart = np
         if (nperod .eq. 3) jstart = np+1
         do 189 j=jstart,n,k
            jm1 = j-np
            jp1 = j+np
            jp2 = jp1+np
            if (jm1 .eq. 0) go to 173
            do 172 i=1,m
               b(i) = q(i,jp1)+q(i,jm1)
  172       continue
            go to 177
  173       if (nperod .eq. 0) go to 175
            do 174 i=1,m
               b(i) = q(i,jp1)
  174       continue
            go to 177
  175       do 176 i=1,m
               b(i) = q(i,jp1)+q(i,n)
  176       continue
  177       call tri (k3,k4,1,m,mm1,ba,bb,bc,b,twocos,d,w)
            do 178 i=1,m
               q(i,j) = q(i,j)+b(i)
  178       continue
            if (jp2 .gt. n) go to 180
            do 179 i=1,m
               b(i) = q(i,j)+q(i,jp2)
  179       continue
            go to 182
  180       do 181 i=1,m
               b(i) = q(i,j)
  181       continue
  182       call tri (k1,k2,1,m,mm1,ba,bb,bc,b,twocos,d,w)
            do 183 i=1,m
               q(i,j) = q(i,j)+b(i)
  183       continue
            if (jp2 .gt. n) go to 185
            do 184 i=1,m
               b(i) = q(i,j)+q(i,jp2)
  184       continue
            go to 187
  185       do 186 i=1,m
               b(i) = q(i,j)
  186       continue
  187       call tri (k3,k4,1,m,mm1,ba,bb,bc,b,twocos,d,w)
            do 188 i=1,m
               q(i,jp1) = q(i,jp1)+b(i)
  188       continue
  189    continue
  190 continue
c
c     regular back substitution for powers of two
c
  191 if (iex2 .eq. 0) go to 202
      np = n2pw
      do 201 kount=1,iex2
         k = np
         np = np/2
         k3 = k-1
         jstart = np
         if (nperod .eq. 3) jstart = 1+np
         do 200 j=jstart,n,k
            jm1 = j-np
            jp1 = j+np
            if (jm1 .ne. 0) go to 194
            if (jp1 .gt. n) go to 200
            if (nperod.eq.1 .or. nperod.eq.2) go to 192
            jm1 = n
            go to 196
  192       do 193 i=1,m
               b(i) = q(i,jp1)
  193       continue
            go to 198
  194       if (jp1 .le. n) go to 196
            do 195 i=1,m
               b(i) = q(i,jm1)
  195       continue
            go to 198
  196       do 197 i=1,m
               b(i) = q(i,jm1)+q(i,jp1)
  197       continue
  198       call tri (np,k3,1,m,mm1,ba,bb,bc,b,twocos,d,w)
            do 199 i=1,m
               q(i,j) = q(i,j)+b(i)
  199       continue
  200    continue
  201 continue
  202 return
      end
      subroutine poinit (nperod,n,iex2,iex3,iex5,n2pw,n2p3p,n2p3p5,ks3,
     1                   ks5,npstop,twocos)
      dimension       twocos(*)
c
c     parameter npstop is not used in this subroutine.
c
c
c     machine dependent constant
c
c     pi=3.1415926535897932384626433832795028841971693993751058209749446
c
      pi = 3.14159265358979
c
c     compute exponents of 2,3, and 5 in n.
c
      np = n
      if (nperod .eq. 1) np = np+1
      if (nperod .eq. 3) np = np-1
      iex2 = 0
  101 k = np/2
      if (2*k .ne. np) go to 102
      iex2 = iex2+1
      np = k
      go to 101
  102 iex3 = 0
  103 k = np/3
      if (3*k .ne. np) go to 104
      iex3 = iex3+1
      np = k
      go to 103
  104 iex5 = 0
  105 k = np/5
      if (5*k .ne. np) go to 106
      iex5 = iex5+1
      np = k
      go to 105
  106 continue
      n2pw = 2**iex2
      n2p3p = n2pw*(3**iex3)
      n2p3p5 = n2p3p*(5**iex5)
c
c     compute necessary values of 2*cos(x)
c
      np = 1
      twocos(1) = 0.
      k = 1
      if (iex2 .eq. 0) go to 110
      l = iex2
      if (iex3+iex5 .ne. 0) go to 107
      if (nperod .eq. 1) l = l-1
      if (l .eq. 0) go to 129
  107 do 109 kount=1,l
         np = 2*np
         do 108 i=1,np
            j = k+i
            twocos(j) = 2.*cos((float(i)-.5)*pi/float(np))
  108    continue
         k = k+np
  109 continue
  110 if (iex3 .eq. 0) go to 114
      l = iex3
      if (iex5 .ne. 0) go to 111
      if (nperod .eq. 1) l = l-1
      if (l .eq. 0) go to 117
  111 do 113 kount=1,l
         np = 3*np
         do 112 i=1,np
            j = k+i
            twocos(j) = 2.*cos((float(i)-.5)*pi/float(np))
  112    continue
         k = k+np
  113 continue
  114 l = iex5
      if (nperod .eq. 1) l = l-1
      if (l .le. 0) go to 117
      do 116 kount=1,l
         np = 5*np
         do 115 i=1,np
            j = k+i
            twocos(j) = 2.*cos((float(i)-.5)*pi/float(np))
  115    continue
         k = k+np
  116 continue
  117 if (nperod.eq.1 .or. nperod.eq.2) go to 121
      if (nperod .eq. 0) go to 119
      npt2 = 2*np
      do 118 i=1,npt2
         j = k+i
         twocos(j) = 2.*cos(float(i)*pi/float(np))
  118 continue
      k = k+npt2
      go to 121
  119 do 120 i=1,np
         j = k+i
         twocos(j) = 2.*cos(2.*float(i)*pi/float(np))
  120 continue
      k = k+np
  121 np = n2p3p5
      if (iex5 .eq. 0) go to 126
      ks5 = k
      do 125 kount=1,iex5
         np = np/5
         npt2 = 2*np
         do 122 i=1,npt2
            j = k+i
            twocos(j) = 2.*cos((float(i)-.5)*pi/float(npt2))
  122    continue
         k = k+npt2
         do 123 i=1,np
            j = k+4*i
            twocos(j-3) = 2.*cos((float(i)-.8)*pi/float(np))
            twocos(j-2) = 2.*cos((float(i)-.6)*pi/float(np))
            twocos(j-1) = 2.*cos((float(i)-.4)*pi/float(np))
            twocos(j) = 2.*cos((float(i)-.2)*pi/float(np))
  123    continue
         k = k+4*np
         do 124 i=1,np
            j = k+2*i
            twocos(j-1) = 2.*cos(float(3*i-2)*pi/float(3*np))
            twocos(j) = 2.*cos(float(3*i-1)*pi/float(3*np))
  124    continue
         k = k+2*np
  125 continue
  126 if (iex3 .eq. 0) go to 129
      ks3 = k
      do 128 kount=1,iex3
         np = np/3
         do 127 i=1,np
            j = k+2*i
            twocos(j-1) = 2.*cos(float(3*i-2)*pi/float(3*np))
            twocos(j) = 2.*cos(float(3*i-1)*pi/float(3*np))
  127    continue
         k = k+2*np
  128 continue
  129 return
      end
      subroutine trid (kstart,kstop,ising,m,mm1,a,b,c,y,twocos,d,w)
      dimension       a(*)       ,b(*)       ,c(*)       ,y(*)       ,
     1                twocos(*)  ,d(*)       ,w(*)
c
c     subroutine to solve tridiagonal systems
c
c
c     parameter w not used in this subroutine.
c
      do 103 k=kstart,kstop
         x = -twocos(k)
         d(1) = c(1)/(b(1)+x)
         y(1) = y(1)/(b(1)+x)
         do 101 i=2,m
            z = b(i)+x-a(i)*d(i-1)
            d(i) = c(i)/z
            y(i) = (y(i)-a(i)*y(i-1))/z
  101    continue
         do 102 ip=1,mm1
            i = m-ip
            y(i) = y(i)-d(i)*y(i+1)
  102    continue
  103 continue
      if (ising .eq. 1) return
      d(1) = c(1)/(b(1)-2.)
      y(1) = y(1)/(b(1)-2.)
      do 104 i=2,mm1
         z = b(i)-2.-a(i)*d(i-1)
         d(i) = c(i)/z
         y(i) = (y(i)-a(i)*y(i-1))/z
  104 continue
      z = b(m)-2.-a(m)*d(m-1)
      x = y(m)-a(m)*y(m-1)
      if (z .ne. 0.) go to 105
      y(m) = 0.
      go to 106
  105 y(m) = x/z
  106 do 107 ip=1,mm1
         i = m-ip
         y(i) = y(i)-d(i)*y(i+1)
  107 continue
      return
      end
      subroutine tridp (kstart,kstop,ising,m,mm1,a,b,c,y,twocos,d,w)
      dimension       a(*)       ,b(*)       ,c(*)       ,y(*)       ,
     1                twocos(*)  ,d(*)       ,w(*)   
c
c     subroutine to solve periodic tridiagonal system
c
      do 103 k=kstart,kstop
         x = -twocos(k)
         d(1) = c(1)/(b(1)+x)
         w(1) = a(1)/(b(1)+x)
         y(1) = y(1)/(b(1)+x)
         bm = b(m)
         z = c(m)
         do 101 i=2,mm1
            den = b(i)+x-a(i)*d(i-1)
            d(i) = c(i)/den
            w(i) = -a(i)*w(i-1)/den
            y(i) = (y(i)-a(i)*y(i-1))/den
            y(m) = y(m)-z*y(i-1)
            bm = bm-z*w(i-1)
            z = -z*d(i-1)
  101    continue
         d(mm1) = d(mm1)+w(mm1)
         z = a(m)+z
         den = bm+x-z*d(mm1)
         y(m) = y(m)-z*y(m-1)
         y(m) = y(m)/den
         y(mm1) = y(mm1)-d(mm1)*y(m)
         do 102 ip=2,mm1
            i = m-ip
            y(i) = y(i)-d(i)*y(i+1)-w(i)*y(m)
  102    continue
  103 continue
      if (ising .eq. 1) return
      d(1) = c(1)/(b(1)-2.)
      w(1) = a(1)/(b(1)-2.)
      y(1) = y(1)/(b(1)-2.)
      bm = b(m)
      z = c(m)
      do 104 i=2,mm1
         den = b(i)-2.-a(i)*d(i-1)
         d(i) = c(i)/den
         w(i) = -a(i)*w(i-1)/den
         y(i) = (y(i)-a(i)*y(i-1))/den
         y(m) = y(m)-z*y(i-1)
         bm = bm-z*w(i-1)
         z = -z*d(i-1)
  104 continue
      d(mm1) = d(mm1)+w(mm1)
      z = a(m)+z
      den = bm-2.-z*d(mm1)
      y(m) = y(m)-z*y(m-1)
      if (den .ne. 0.) go to 105
      y(m) = 0.
      go to 106
  105 y(m) = y(m)/den
  106 y(mm1) = y(mm1)-d(mm1)*y(m)
      do 107 ip=2,mm1
         i = m-ip
         y(i) = y(i)-d(i)*y(i+1)-w(i)*y(m)
  107 continue
      return
      end
      function apxeps (dum)
c          return machine acc.
      apxeps = 1.00  e -4
      do 20 i = 1,100
      s = apxeps * 0.01
      s = s + 1.
      s = s - 1.
c     print 10,i,s,apxeps
      if (s.le.dum) return
   20 apxeps = apxeps *.1
      return
   10 format (10x,'i = ',i4,' s = ',f20.18,' eps = ',f20.18)
      end

