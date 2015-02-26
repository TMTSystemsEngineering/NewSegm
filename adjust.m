function f = adjust(gamma,supn,Rnom,n_segments);
        R=Rnom*gamma;
        t=[90 210 330];
        SSASup_PSA=[R*cosd(t); R*sind(t)];   % each column contains coordinates of one SSA support point in the PSA system
        % calculate residual adjustments of AAP
        for i=1:n_segments
            % calculate cell node to SSA support distances (required adjustment of AAP)
            a=supn(1:2,:,i)-SSASup_PSA;
            aa(:,i)=sqrt( a(1,:).^2 + a(2,:).^2 );
        end;
        f = max(max(aa));