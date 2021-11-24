function  [Ln,labelstot] = UGM_Decode_ICM_3D_v1(labels,dist,dist_co,neighbor,griddim,nodePot,level,k_down,useTop,previous,x_tot,y_tot,z_tot,labelstot,grid_space,spc,quant,smooth_co_tot,Tcoe_n_tot,top_co_tot)
num_neighbor=6;
[junk Ln] = min(nodePot,[],2);
nneighbor=size(neighbor,2);
sample_space_x=quant(level,1);
sample_space_y=quant(level,2);
sample_space_z=quant(level,3);
top_co=top_co_tot(level);
Tcoe_n=Tcoe_n_tot(level);
r=griddim(1);
c=griddim(2);
s=griddim(3);
x_h=fix(x_tot/2)+1;
y_h=fix(y_tot/2)+1;
z_h=fix(z_tot/2)+1;

% itsmooth=size(neighbor,2);

% top_co=spc(1)*spc(2)*spc(3)/(grid_space(level,1))^3;
% top_co=top_co_tot(level);
% Tcoe_n=Tcoe_n_tot(level);

cur_grid_space_x=grid_space(level,1);
cur_grid_space_y=grid_space(level,2);
cur_grid_space_z=grid_space(level,3);
cur_pix_spc=[cur_grid_space_x,cur_grid_space_y,cur_grid_space_z];
coe_x=10/cur_pix_spc(1);
coe_y=10/cur_pix_spc(2);
coe_z=10/cur_pix_spc(3);
% spc_tot=ones(itsmooth,3).*spc;

smooth_co=smooth_co_tot(level);
[x_index_tot,y_index_tot,z_index_tot]=Label_Coordinate_tot(x_tot,y_tot,z_tot);

z_index_tot=-z_index_tot;

x_index_tot_spc=x_index_tot.*spc(1);
y_index_tot_spc=y_index_tot.*spc(2);
z_index_tot_spc=z_index_tot.*spc(3);


labelindex=1:x_tot*y_tot*z_tot;
labelindex=reshape(labelindex,[x_tot,y_tot,z_tot]);



% This label index the data term
labelnum=zeros(r*c*s+1,1)+floor(labels.nlabels/2)+1;
labelstot_1=zeros(r*c*s+1,1)+floor(x_tot*y_tot*z_tot/2)+1;
labelstot_1(1:r*c*s)=labelstot;
labelstot=labelstot_1;
  
% totenergyend=zeros(5000,1);
maxIter=50;
% niter=fix(maxIter*k_down^(-level+7));
% randnewlabel=unidrnd(labels.nlabels,[r*c*s*niter,1]);
[x_index,y_index,z_index]=Label_Coordinate(labels);
x_index=x_index*sample_space_x;
y_index=y_index*sample_space_y;
z_index=z_index*sample_space_z;



[nNodes,maxStates] = size(nodePot);

nStates = maxStates;
y=labelnum;
done = 0;
iter=1;
while ~done&&(iter<maxIter)
    done = 1;
	y2 = y;
    for n = 1:nNodes
        
        x_mov_or_add=x_index+previous(n,1)+x_h;
        y_mov_or_add=y_index+previous(n,2)+y_h;
        z_mov_or_add=z_index+previous(n,3)+z_h;
        totlabel=zeros(1,nStates);
        for idx=1:nStates 
            totlabel(idx)=labelindex(x_mov_or_add(idx),y_mov_or_add(idx),z_mov_or_add(idx));
        end
    
        % Compute Node Potential
        pot = nodePot(n,1:nStates);

        % Find Neighbors
        edges=neighbor(n,:);
        if num_neighbor==6
            up_point=edges(1);
            down_point=edges(6);
            left_point=edges(2);
            right_point=edges(5);
            anter_point=edges(3);
            poster_point=edges(4);
        else
            if num_neighbor==18
                up_point=edges(3);
                down_point=edges(16);
                left_point=edges(7);
                right_point=edges(12);
                anter_point=edges(9);
                poster_point=edges(10);
            end
        end
                
        

            smooth_cur=zeros(nStates,3);
            x_index_nstate=x_index_tot_spc(totlabel);
            y_index_nstate=y_index_tot_spc(totlabel);
            z_index_nstate=z_index_tot_spc(totlabel);
            smooth_cur(:,1)=x_index_nstate;
            smooth_cur(:,2)=y_index_nstate;
            smooth_cur(:,3)=z_index_nstate;
            eptot=zeros(1,nStates);

            for e=1:nneighbor
                if isempty(dist)
                ecurr=labelstot(edges(e));
                ecurr_x=x_index_tot_spc(ecurr);
                ecurr_y=y_index_tot_spc(ecurr);
                ecurr_z=z_index_tot_spc(ecurr);
                smooth_e=repmat([ecurr_x,ecurr_y,ecurr_z],nStates,1);
                
%                 ep=sum(abs(smooth_cur-smooth_e),2);
                ep=(sum((smooth_cur-smooth_e).^2,2));
                ep(find(ep>=900))=900;
%                 ep=ep./dist_co/3;
                ep=ep';
                eptot=eptot+ep;

%                 ep=zeros(1,nStates);
%                 pot=pot+ep;
%                 orb_tot(orb_tot>=80)=80;
%                 newb_tot=sqrt(sum((smooth_new-smooth_neighbor).^2,2));
%                 newb_tot(newb_tot>=80)=80;
%                 orbenergy=sum(orb_tot)/(dist_co)*smooth_co;
%                 newbenergy=sum(newb_tot)/(dist_co)*smooth_co; 
                else
                    ecurr=labelstot(edges(e));
                    ep=dist(ecurr,totlabel);
                    ep=ep.*smooth_co;
                    pot = pot +ep;
                end
            end

pot=pot+sqrt(eptot)./dist_co*3*smooth_co;

        % Multiply Edge Potentials

        etp=zeros(1,nStates);

        if useTop==1
            for tpidx=1:nStates
                ortotlabel=totlabel(tpidx);
                tporpx=x_index_tot_spc(ortotlabel);
                tporpy=y_index_tot_spc(ortotlabel);
                tporpz=z_index_tot_spc(ortotlabel);
    
                or_p_up_x=x_index_tot_spc(labelstot(up_point));
                or_p_up_y=y_index_tot_spc(labelstot(up_point));
                or_p_up_z=z_index_tot_spc(labelstot(up_point));
                or_p_up_z=(or_p_up_z+cur_grid_space_z*spc(3));
             

                or_p_down_x=x_index_tot_spc(labelstot(down_point));
                or_p_down_y=y_index_tot_spc(labelstot(down_point));
                or_p_down_z=z_index_tot_spc(labelstot(down_point));
                or_p_down_z=(or_p_down_z-cur_grid_space_z*spc(3));
        

                or_p_left_x=x_index_tot_spc(labelstot(left_point));
                or_p_left_y=y_index_tot_spc(labelstot(left_point));
                or_p_left_z=z_index_tot_spc(labelstot(left_point));
                or_p_left_y=(or_p_left_y-cur_grid_space_y*spc(2));
        
                or_p_right_x=x_index_tot_spc(labelstot(right_point));
                or_p_right_y=y_index_tot_spc(labelstot(right_point));
                or_p_right_z=z_index_tot_spc(labelstot(right_point));
                or_p_right_y=or_p_right_y+cur_grid_space_y*spc(2);
        

                or_p_anter_x=x_index_tot_spc(labelstot(anter_point));
                or_p_anter_y=y_index_tot_spc(labelstot(anter_point));
                or_p_anter_z=z_index_tot_spc(labelstot(anter_point));
                or_p_anter_x=or_p_anter_x-cur_grid_space_x*spc(1);
        
         

                or_p_poster_x=x_index_tot_spc(labelstot(poster_point));
                or_p_poster_y=y_index_tot_spc(labelstot(poster_point));
                or_p_poster_z=z_index_tot_spc(labelstot(poster_point)); 
                or_p_poster_x=or_p_poster_x+cur_grid_space_x*spc(1);
                
                

                 %Jaco1 point 1 4 5
%                  Jaco_1=(orpx-or_p_anter_x)*(or_p_right_y-orpy)*(or_p_up_z-orpz)+(orpy-or_p_anter_y)*(or_p_right_z-orpz)*(or_p_up_x-orpx)+(or_p_right_x-orpx)*(or_p_up_y-orpy)*(orpz-or_p_anter_z)-(orpz-or_p_anter_z)*(or_p_right_y-orpy)*(or_p_up_x-orpx)-(orpy-or_p_anter_y)*(or_p_right_x-orpx)*(or_p_up_z-orpz)-(or_p_right_z-orpz)*(or_p_up_y-orpy)*(orpx-or_p_anter_x);
                Jaco_1=(tporpx-or_p_anter_x)*(or_p_right_y-tporpy)*(or_p_up_z-tporpz)+(tporpy-or_p_anter_y)*(or_p_right_z-tporpz)*(or_p_up_x-tporpx)+(or_p_right_x-tporpx)*(or_p_up_y-tporpy)*(tporpz-or_p_anter_z)-(tporpz-or_p_anter_z)*(or_p_right_y-tporpy)*(or_p_up_x-tporpx)-(tporpy-or_p_anter_y)*(or_p_right_x-tporpx)*(or_p_up_z-tporpz)-(or_p_right_z-tporpz)*(or_p_up_y-tporpy)*(tporpx-or_p_anter_x);
                Jaco_1=Jaco_1*top_co;
                if Jaco_1>=0
                    orthird_1=0;
                else
                    orthird_1=log(-Jaco_1+1)*Tcoe_n;
                end
        
        
                %Jaco2 point 1 3 5
%                  Jaco_2=(orpx-or_p_anter_x)*(orpy-or_p_left_y)*(or_p_up_z-orpz)+(orpy-or_p_anter_y)*(orpz-or_p_left_z)*(or_p_up_x-orpx)+(orpx-or_p_left_x)*(or_p_up_y-orpy)*(orpz-or_p_anter_z)-(orpz-or_p_anter_z)*(orpy-or_p_left_y)*(or_p_up_x-orpx)-(orpy-or_p_anter_y)*(orpx-or_p_left_x)*(or_p_up_z-orpz)-(orpz-or_p_left_z)*(or_p_up_y-orpy)*(orpx-or_p_anter_x);
                Jaco_2=(tporpx-or_p_anter_x)*(tporpy-or_p_left_y)*(or_p_up_z-tporpz)+(tporpy-or_p_anter_y)*(tporpz-or_p_left_z)*(or_p_up_x-tporpx)+(tporpx-or_p_left_x)*(or_p_up_y-tporpy)*(tporpz-or_p_anter_z)-(tporpz-or_p_anter_z)*(tporpy-or_p_left_y)*(or_p_up_x-tporpx)-(tporpy-or_p_anter_y)*(tporpx-or_p_left_x)*(or_p_up_z-tporpz)-(tporpz-or_p_left_z)*(or_p_up_y-tporpy)*(tporpx-or_p_anter_x);
                Jaco_2=Jaco_2*top_co;
                if Jaco_2>=0
                   orthird_2=0;
                else
                   orthird_2=log(-Jaco_2+1)*Tcoe_n;
                end
        
        
                %Jaco3 point 2 4 5
%                 Jaco_3=(orpx-or_p_anter_x)*(or_p_right_y-orpy)*(orpz-or_p_down_z)+(orpy-or_p_anter_y)*(or_p_right_z-orpz)*(orpx-or_p_down_x)+(or_p_right_x-orpx)*(orpy-or_p_down_y)*(orpz-or_p_anter_z)-(orpz-or_p_anter_z)*(or_p_right_y-orpy)*(orpx-or_p_down_x)-(orpy-or_p_anter_y)*(or_p_right_x-orpx)*(orpz-or_p_down_z)-(or_p_right_z-orpz)*(orpy-or_p_down_y)*(orpx-or_p_anter_x); 
                Jaco_3=(tporpx-or_p_anter_x)*(or_p_right_y-tporpy)*(tporpz-or_p_down_z)+(tporpy-or_p_anter_y)*(or_p_right_z-tporpz)*(tporpx-or_p_down_x)+(or_p_right_x-tporpx)*(tporpy-or_p_down_y)*(tporpz-or_p_anter_z)-(tporpz-or_p_anter_z)*(or_p_right_y-tporpy)*(tporpx-or_p_down_x)-(tporpy-or_p_anter_y)*(or_p_right_x-tporpx)*(tporpz-or_p_down_z)-(or_p_right_z-tporpz)*(tporpy-or_p_down_y)*(tporpx-or_p_anter_x);
                Jaco_3=Jaco_3*top_co;
        
                if Jaco_3>=0
                    orthird_3=0;
                else


                    orthird_3=log(-Jaco_3+1)*Tcoe_n;
                end
        
        
                %Jaco4 point 2 3 5
%                                Jaco_4=(orpx-or_p_anter_x)*(orpy-or_p_left_y)*(orpz-or_p_down_z)+(orpy-or_p_anter_y)*(orpz-or_p_left_z)*(orpx-or_p_down_x)+(orpx-or_p_left_x)*(orpy-or_p_down_y)*(orpz-or_p_anter_z)-(orpz-or_p_anter_z)*(orpy-or_p_left_y)*(orpx-or_p_down_x)-(orpy-or_p_anter_y)*(orpx-or_p_left_x)*(orpz-or_p_down_z)-(orpz-or_p_left_z)*(orpy-or_p_down_y)*(orpx-or_p_anter_x);
                Jaco_4=(tporpx-or_p_anter_x)*(tporpy-or_p_left_y)*(tporpz-or_p_down_z)+(tporpy-or_p_anter_y)*(tporpz-or_p_left_z)*(tporpx-or_p_down_x)+(tporpx-or_p_left_x)*(tporpy-or_p_down_y)*(tporpz-or_p_anter_z)-(tporpz-or_p_anter_z)*(tporpy-or_p_left_y)*(tporpx-or_p_down_x)-(tporpy-or_p_anter_y)*(tporpx-or_p_left_x)*(tporpz-or_p_down_z)-(tporpz-or_p_left_z)*(tporpy-or_p_down_y)*(tporpx-or_p_anter_x);
                Jaco_4=Jaco_4*top_co; 
        
                if Jaco_4>=0
                    orthird_4=0;
                else
                    orthird_4=log(-Jaco_4+1)*Tcoe_n;
                end
        
        
                %Jaco5 point 1 4 6 fff
                %Jaco_5=(p_right_x-orpx)*(p_poster_y-orpy)*(p_up_z-orpz)+(p_poster_x-orpx)*(p_up_y-orpy)*(p_right_z-orpz)+(p_right_y-orpy)*(p_poster_z-orpz)*(p_up_x-orpx)-(p_right_z-orpz)*(p_poster_y-orpy)*(p_up_x-orpx)-(p_right_y-orpy)*(p_poster_x-orpx)*(p_up_z-orpz)-(p_poster_z-orpz)*(p_up_y-orpy)*(p_right_x-orpx);
%                              Jaco_5=(or_p_poster_x-orpx)*(or_p_right_y-orpy)*(or_p_up_z-orpz)+(or_p_poster_y-orpy)*(or_p_right_z-orpz)*(or_p_up_x-orpx)+(or_p_right_x-orpx)*(or_p_up_y-orpy)*(or_p_poster_z-orpz)-(or_p_poster_z-orpz)*(or_p_right_y-orpy)*(or_p_up_x-orpx)-(or_p_poster_y-orpy)*(or_p_right_x-orpx)*(or_p_up_z-orpz)-(or_p_right_z-orpz)*(or_p_up_y-orpy)*(or_p_poster_x-orpx);
                Jaco_5=(or_p_poster_x-tporpx)*(or_p_right_y-tporpy)*(or_p_up_z-tporpz)+(or_p_poster_y-tporpy)*(or_p_right_z-tporpz)*(or_p_up_x-tporpx)+(or_p_right_x-tporpx)*(or_p_up_y-tporpy)*(or_p_poster_z-tporpz)-(or_p_poster_z-tporpz)*(or_p_right_y-tporpy)*(or_p_up_x-tporpx)-(or_p_poster_y-tporpy)*(or_p_right_x-tporpx)*(or_p_up_z-tporpz)-(or_p_right_z-tporpz)*(or_p_up_y-tporpy)*(or_p_poster_x-tporpx);
                Jaco_5=Jaco_5*top_co;
                if Jaco_5>=0
                   orthird_5=0;
                else

                   orthird_5=log(-Jaco_5+1)*Tcoe_n;
                end
        
        
        
                 %Jaco6 point 2 3 6
                 
%                   Jaco_6=(or_p_poster_x-orpx)*(orpy-or_p_left_y)*(orpz-or_p_down_z)+(or_p_poster_y-orpy)*(orpz-or_p_left_z)*(orpx-or_p_down_x)+(orpx-or_p_left_x)*(orpy-or_p_down_y)*(or_p_poster_z-orpz)-(or_p_poster_z-orpz)*(orpy-or_p_left_y)*(orpx-or_p_down_x)-(or_p_poster_y-orpy)*(orpx-or_p_left_x)*(orpz-or_p_down_z)-(orpz-or_p_left_z)*(orpy-or_p_down_y)*(or_p_poster_x-orpx);
                  %Jaco_6=(orpx-p_left_x)*(p_poster_y-orpy)*(orpz-p_down_z)+(p_poster_x-orpx)*(orpy-p_down_y)*(orpz-p_left_z)+(orpy-p_left_y)*(p_poster_z-orpz)*(orpx-p_down_x)-(orpz-p_left_z)*(p_poster_y-orpy)*(orpx-p_down_x)-(orpy-p_left_y)*(p_poster_x-orpx)*(orpz-p_down_z)-(p_poster_z-orpz)*(orpy-p_down_y)*(orpx-p_left_x);
                 Jaco_6=(or_p_poster_x-tporpx)*(tporpy-or_p_left_y)*(tporpz-or_p_down_z)+(or_p_poster_y-tporpy)*(tporpz-or_p_left_z)*(tporpx-or_p_down_x)+(tporpx-or_p_left_x)*(tporpy-or_p_down_y)*(or_p_poster_z-tporpz)-(or_p_poster_z-tporpz)*(tporpy-or_p_left_y)*(tporpx-or_p_down_x)-(or_p_poster_y-tporpy)*(tporpx-or_p_left_x)*(tporpz-or_p_down_z)-(tporpz-or_p_left_z)*(tporpy-or_p_down_y)*(or_p_poster_x-tporpx);
                 Jaco_6=Jaco_6*top_co;
                 if Jaco_6>=0
                     orthird_6=0;
                 else

                    orthird_6=log(-Jaco_6+1)*Tcoe_n;
                 end
        
        
        
               %Jaco7 point 2 4 6
               %Jaco_7=(p_right_x-orpx)*(p_poster_y-orpy)*(orpz-p_down_z)+(p_poster_x-orpx)*(orpy-p_down_y)*(p_right_z-orpz)+(p_right_y-orpy)*(p_poster_z-orpz)*(orpx-p_down_x)-(p_right_z-orpz)*(p_poster_y-orpy)*(orpx-p_down_x)-(p_right_y-orpy)*(p_poster_x-orpx)*(orpz-p_down_z)-(p_poster_z-orpz)*(orpy-p_down_y)*(-p_right_x-orpx);
%                      Jaco_7=(or_p_poster_x-orpx)*(or_p_right_y-orpy)*(orpz-or_p_down_z)+(or_p_poster_y-orpy)*(or_p_right_z-orpz)*(orpx-or_p_down_x)+(or_p_right_x-orpx)*(orpy-or_p_down_y)*(or_p_poster_z-orpz)-(or_p_poster_z-orpz)*(or_p_right_y-orpy)*(orpx-or_p_down_x)-(or_p_poster_y-orpy)*(or_p_right_x-orpx)*(orpz-or_p_down_z)-(or_p_right_z-orpz)*(orpy-or_p_down_y)*(or_p_poster_x-orpx);
               Jaco_7=(or_p_poster_x-tporpx)*(or_p_right_y-tporpy)*(tporpz-or_p_down_z)+(or_p_poster_y-tporpy)*(or_p_right_z-tporpz)*(tporpx-or_p_down_x)+(or_p_right_x-tporpx)*(tporpy-or_p_down_y)*(or_p_poster_z-tporpz)-(or_p_poster_z-tporpz)*(or_p_right_y-tporpy)*(tporpx-or_p_down_x)-(or_p_poster_y-tporpy)*(or_p_right_x-tporpx)*(tporpz-or_p_down_z)-(or_p_right_z-tporpz)*(tporpy-or_p_down_y)*(or_p_poster_x-tporpx);
               Jaco_7=Jaco_7*top_co;
               if Jaco_7>=0
                   orthird_7=0;
               else

                   orthird_7=log(-Jaco_7+1)*Tcoe_n;
               end
        
               %Jaco8 point 1 3 6
              
               %Jaco_8=(orpx-p_left_x)*(p_poster_y-orpy)*(p_up_z-orpz)+(p_poster_x-orpx)*(p_up_y-orpy)*(orpz-p_left_z)+(orpy-p_left_y)*(p_poster_z-orpz)*(p_up_x-orpx)-(orpz-p_left_z)*(p_poster_y-orpy)*(p_up_x-orpx)-(orpy-p_left_y)*(p_poster_x-orpx)*(p_up_z-orpz)-(p_poster_z-orpz)*(p_up_y-orpy)*(orpx-p_left_x);
%                          Jaco_8=(or_p_poster_x-orpx)*(orpy-or_p_left_y)*(or_p_up_z-orpz)+(or_p_poster_y-orpy)*(orpz-or_p_left_z)*(or_p_up_x-orpx)+(orpx-or_p_left_x)*(or_p_up_y-orpy)*(or_p_poster_z-orpz)-(or_p_poster_z-orpz)*(orpy-or_p_left_y)*(or_p_up_x-orpx)-(or_p_poster_y-orpy)*(orpx-or_p_left_x)*(or_p_up_z-orpz)-(orpz-or_p_left_z)*(or_p_up_y-orpy)*(or_p_poster_x-orpx);
               Jaco_8=(or_p_poster_x-tporpx)*(tporpy-or_p_left_y)*(or_p_up_z-tporpz)+(or_p_poster_y-tporpy)*(tporpz-or_p_left_z)*(or_p_up_x-tporpx)+(tporpx-or_p_left_x)*(or_p_up_y-tporpy)*(or_p_poster_z-tporpz)-(or_p_poster_z-tporpz)*(tporpy-or_p_left_y)*(or_p_up_x-tporpx)-(or_p_poster_y-tporpy)*(tporpx-or_p_left_x)*(or_p_up_z-tporpz)-(tporpz-or_p_left_z)*(or_p_up_y-tporpy)*(or_p_poster_x-tporpx);
               Jaco_8=Jaco_8*top_co;
        
        
                if Jaco_8>=0
                    orthird_8=0;
                else
%                     orthird_8=1000;
                    orthird_8=log(-Jaco_8+1)*Tcoe_n;
                end
                etp(tpidx)=orthird_1+orthird_2+orthird_3+orthird_4+orthird_5+orthird_6+orthird_7+orthird_8;
        

            end
        end
        

        pot=pot+etp;
        % Assign to Maximum State
        [junk newY] = min(pot);
        if newY ~= y(n)
            y(n) = newY;
            labelstot(n)=totlabel(newY);
%             done = 0;
        end
    end
   iter=iter+1;
   changes=sum(y2~=y);
   if changes~=0
       done=0;
   else
       done=1;
   end
%    fprintf('logPot = %f, changes = %d\n',UGM_ConfigurationPotential_correct_TP(y,nodePot,edgePot,edgeEnds,TP_energy,C_tp),sum(y2~=y));
   
   fprintf('changes = %d, iter = %d\n',sum(y2~=y), iter);
%             fprintf('changes = %d\n',sum(y2~=y));
 end


Ln=y(1:nNodes);
end
