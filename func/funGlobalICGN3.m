% =========================================================
% function GlobalICGN to solve FE iterations for global DVC
% ---------------------------------------------------------
%   INPUT:
%       DVC mesh, DVC image pair,
%       displacement initial guess, regularizer coefficient
%
%   OUTPUT:
%       U: Solved displacement field;
%       normOfW: FE-based global DVC iteration update norm;
%       timeICGN: Time cost for each FE-based global DVC iteration;
%
% Author: Jin Yang, jyang526@wisc.edu or aldicdvc@gmail.com
% Date: 2020.10
% =========================================================

function [U,normOfW,timeICGN] = funGlobalICGN3(DVCmesh,Df,Img1,Img2,U,alpha,tol,maxIter)

coordinatesFEM = DVCmesh.coordinatesFEM; % FE-mesh coordinates
elementsFEM = DVCmesh.elementsFEM; % FE-mesh element
DIM = 3; % Problem dimension
NodesPerEle = 8; % Using cubic elements
FEMSize = DIM*size(coordinatesFEM,1); % FE-system size
winsize = (coordinatesFEM(2,1)-coordinatesFEM(1,1))*ones(1,3); % Finite element size

DfDx = Df.DfDx; DfDy = Df.DfDy; DfDz = Df.DfDz; % ROI image parameters
DfAxis = Df.DfAxis; DfDxStartx = DfAxis(1); DfDxStarty = DfAxis(3); DfDxStartz = DfAxis(5);
try maxIter = maxIter; catch maxIter = 100; end % set max iteration number as 100 by default


%% %%%%%%%%%%%%%%%% Start FE ICGN iteration %%%%%%%%%%%%%%%
for stepwithinwhile = 1:maxIter
    
    tic;
    
    % ====== Initialize stiffness matrix at the first iteration ======
    if (stepwithinwhile==1)
        disp(['--- Global IC-GN iterations ---']);
        INDEXAI = []; INDEXAJ = []; INDEXAVAL = []; INDEXAREG = []; % A = sparse(FEMSize,FEMSize);
    end
    INDEXBI = []; INDEXBVAL = []; % clear b; b = sparse(FEMsize,1); % Update external force vector
    
    
    
    % ------ Define ksi, eta and zeta list ------
    ksiList = -1:2/winsize(1):1; etaList = -1:2/winsize(2):1; zetaList = -1:2/winsize(3):1;
    [ksiMat,etaMat,zetaMat] = ndgrid(ksiList,etaList,zetaList);
    
    NMat = cell(8,1);
    NMat{1} = 1/8*(1-ksiMat).*(1-etaMat).*(1-zetaMat); NMat{2} = 1/8*(1+ksiMat).*(1-etaMat).*(1-zetaMat);
    NMat{3} = 1/8*(1+ksiMat).*(1+etaMat).*(1-zetaMat); NMat{4} = 1/8*(1-ksiMat).*(1+etaMat).*(1-zetaMat);
    NMat{5} = 1/8*(1-ksiMat).*(1-etaMat).*(1+zetaMat); NMat{6} = 1/8*(1+ksiMat).*(1-etaMat).*(1+zetaMat);
    NMat{7} = 1/8*(1+ksiMat).*(1+etaMat).*(1+zetaMat); NMat{8} = 1/8*(1-ksiMat).*(1+etaMat).*(1+zetaMat);
    
    hbar = waitbar(0,['Global ICGN iteartion step: ',num2str(stepwithinwhile)]);
    % ============= Each element, assemble stiffness matrix ============
    for indEle = 1 : size(elementsFEM,1) % indEle is the element index
        
        waitbar(indEle/size(elementsFEM,1));
        
        tempA = zeros(DIM*NodesPerEle,DIM*NodesPerEle); tempb = tempA(:,1);
        
        % ------ Find corner pts ------
        pt1xyz = coordinatesFEM(elementsFEM(indEle,1),:); pt2xyz = coordinatesFEM(elementsFEM(indEle,2),:);
        pt3xyz = coordinatesFEM(elementsFEM(indEle,3),:); pt4xyz = coordinatesFEM(elementsFEM(indEle,4),:);
        pt5xyz = coordinatesFEM(elementsFEM(indEle,5),:); pt6xyz = coordinatesFEM(elementsFEM(indEle,6),:);
        pt7xyz = coordinatesFEM(elementsFEM(indEle,7),:); pt8xyz = coordinatesFEM(elementsFEM(indEle,8),:);
        pt1x = pt1xyz(1); pt1y = pt1xyz(2); pt1z = pt1xyz(3); pt2x = pt2xyz(1); pt2y = pt2xyz(2); pt2z = pt2xyz(3);
        pt3x = pt3xyz(1); pt3y = pt3xyz(2); pt3z = pt3xyz(3); pt4x = pt4xyz(1); pt4y = pt4xyz(2); pt4z = pt4xyz(3);
        pt5x = pt5xyz(1); pt5y = pt5xyz(2); pt5z = pt5xyz(3); pt6x = pt6xyz(1); pt6y = pt6xyz(2); pt6z = pt6xyz(3);
        pt7x = pt7xyz(1); pt7y = pt7xyz(2); pt7z = pt7xyz(3); pt8x = pt8xyz(1); pt8y = pt8xyz(2); pt8z = pt8xyz(3);
        
        % ------ Calculate ksi and eta ------
        % lMatrix = [ pt1x*pt1y*pt1z, pt1x*pt1y, pt1y*pt1z, pt1z*pt1x, pt1x, pt1y, pt1z, 1;
        %             pt2x*pt2y*pt2z, pt2x*pt2y, pt2y*pt2z, pt2z*pt2x, pt2x, pt2y, pt2z, 1;
        %             pt3x*pt3y*pt3z, pt3x*pt3y, pt3y*pt3z, pt3z*pt3x, pt3x, pt3y, pt3z, 1;
        %             pt4x*pt4y*pt4z, pt4x*pt4y, pt4y*pt4z, pt4z*pt4x, pt4x, pt4y, pt4z, 1;
        %             pt5x*pt5y*pt5z, pt5x*pt5y, pt5y*pt5z, pt5z*pt5x, pt5x, pt5y, pt5z, 1;
        %             pt6x*pt6y*pt6z, pt6x*pt6y, pt6y*pt6z, pt6z*pt6x, pt6x, pt6y, pt6z, 1;
        %             pt7x*pt7y*pt7z, pt7x*pt7y, pt7y*pt7z, pt7z*pt7x, pt7x, pt7y, pt7z, 1;
        %             pt8x*pt8y*pt8z, pt8x*pt8y, pt8y*pt8z, pt8z*pt8x, pt8x, pt8y, pt8z, 1];
        
        % ------ Find linear interpolation coefficients ------
        % lb = [-1;1;1;-1;-1;1;1;-1]; lCoeff = linsolve(lMatrix,lb);
        % mb = [-1;-1;1;1;-1;-1;1;1]; mCoeff = linsolve(lMatrix,mb);
        % nb = [-1;-1;-1;-1;1;1;1;1]; nCoeff = linsolve(lMatrix,nb);
        
        % ------ Find element nodal indices ------
        tp = ones(1,DIM);
        tempIndexU = 3*elementsFEM(indEle,[tp,2*tp,3*tp,4*tp,5*tp,6*tp,7*tp,8*tp]);
        tempIndexU(1:3:end) = tempIndexU(1:3:end)-2;
        tempIndexU(2:3:end) = tempIndexU(2:3:end)-1; % size of tempIndexU: 1*24
        
        
        [ptOfxAll,ptOfyAll,ptOfzAll] = ndgrid(pt1x:pt7x,pt1y:pt7y,pt1z:pt7z); % To compute at each pixels
        
        % U1Mat = U(tempIndexU(3*1-2))*ones(winsize+ones(1,3)); U2Mat = U(tempIndexU(3*2-2))*ones(winsize+ones(1,3));
        % U3Mat = U(tempIndexU(3*3-2))*ones(winsize+ones(1,3)); U4Mat = U(tempIndexU(3*4-2))*ones(winsize+ones(1,3));
        % U5Mat = U(tempIndexU(3*5-2))*ones(winsize+ones(1,3)); U6Mat = U(tempIndexU(3*6-2))*ones(winsize+ones(1,3));
        % U7Mat = U(tempIndexU(3*7-2))*ones(winsize+ones(1,3)); U8Mat = U(tempIndexU(3*8-2))*ones(winsize+ones(1,3));
        tempUMat = zeros(winsize+ones(1,3)); tempVMat = tempUMat; tempWMat = tempUMat;
        for tempk = 1:NodesPerEle
            tempUMat = tempUMat + (U(tempIndexU(3*tempk-2))*ones(winsize+ones(1,3))).*NMat{tempk};
            tempVMat = tempVMat + (U(tempIndexU(3*tempk-1))*ones(winsize+ones(1,3))).*NMat{tempk};
            tempWMat = tempWMat + (U(tempIndexU(3*tempk-0))*ones(winsize+ones(1,3))).*NMat{tempk};
        end
        %ptOfxMat = ptOfxAll + tempUMat; ptOfyMat = ptOfyAll + tempVMat; ptOfzMat = ptOfzAll + tempWMat;
        
        tempg = ba_interp3(Img2, ptOfyAll+tempVMat, ptOfxAll+tempUMat, ptOfzAll+tempWMat, 'cubic'); % Deformed g(x+u)
        
        
        ptOfxAll = ptOfxAll(:); ptOfyAll = ptOfyAll(:); ptOfzAll = ptOfzAll(:); % Garantuee ptOfxAll, pyOfyAll, and ptOfzAll are column vectors
        for tempjj = 1:length(ptOfxAll) % Write into one for-loop instead of three for-loops
            % for ptOfx = pt1x:pt7x
            %    for ptOfy = pt1y:pt7y
            %        for ptOfz = pt1z:pt7z
            
            ptOfx = ptOfxAll(tempjj); ptOfy = ptOfyAll(tempjj); ptOfz = ptOfzAll(tempjj);
            
            % ------ Calculate ksi, eta and zeta ------
            ksi = ksiMat(tempjj); eta = etaMat(tempjj); zeta = zetaMat(tempjj);
            %ksi = [ptOfx*ptOfy*ptOfz, ptOfx*ptOfy, ptOfy*ptOfz, ptOfz*ptOfx, ptOfx, ptOfy, ptOfz, 1]*lCoeff;
            %eta = [ptOfx*ptOfy*ptOfz, ptOfx*ptOfy, ptOfy*ptOfz, ptOfz*ptOfx, ptOfx, ptOfy, ptOfz, 1]*mCoeff;
            %zeta = [ptOfx*ptOfy*ptOfz, ptOfx*ptOfy, ptOfy*ptOfz, ptOfz*ptOfx, ptOfx, ptOfy, ptOfz, 1]*nCoeff;
            
            % ------ Calculate N matrix ------
            N1 = NMat{1}(tempjj); N2 = NMat{2}(tempjj); N3 = NMat{3}(tempjj); N4 = NMat{4}(tempjj);
            N5 = NMat{5}(tempjj); N6 = NMat{6}(tempjj); N7 = NMat{7}(tempjj); N8 = NMat{8}(tempjj);
            %N1 = 1/8*(1-ksi)*(1-eta)*(1-zeta); N2 = 1/8*(1+ksi)*(1-eta)*(1-zeta);
            %N3 = 1/8*(1+ksi)*(1+eta)*(1-zeta); N4 = 1/8*(1-ksi)*(1+eta)*(1-zeta);
            %N5 = 1/8*(1-ksi)*(1-eta)*(1+zeta); N6 = 1/8*(1+ksi)*(1-eta)*(1+zeta);
            %N7 = 1/8*(1+ksi)*(1+eta)*(1+zeta); N8 = 1/8*(1-ksi)*(1+eta)*(1+zeta);
            
            % ------ Generate [N] shape function matrix ------
            tpN1 = diag([N1,N1,N1]); tpN2 = diag([N2,N2,N2]); tpN3 = diag([N3,N3,N3]); tpN4 = diag([N4,N4,N4]);
            tpN5 = diag([N5,N5,N5]); tpN6 = diag([N6,N6,N6]); tpN7 = diag([N7,N7,N7]); tpN8 = diag([N8,N8,N8]);
            NOrig = [tpN1,tpN2,tpN3,tpN4,tpN5,tpN6,tpN7,tpN8]; N = NOrig;
            NDiag = diag([N1*ones(1,DIM),N2*ones(1,DIM),N3*ones(1,DIM),N4*ones(1,DIM),N5*ones(1,DIM),N6*ones(1,DIM),N7*ones(1,DIM),N8*ones(1,DIM)]);
            
            % ------ Build J matrix ------
            % Comment: I didn't change Jacobian matrix J when enriched
            % functions are added.
            J = [funDN1Dksi(ksi,eta,zeta),funDN2Dksi(ksi,eta,zeta),funDN3Dksi(ksi,eta,zeta),funDN4Dksi(ksi,eta,zeta), ...
                funDN5Dksi(ksi,eta,zeta),funDN6Dksi(ksi,eta,zeta),funDN7Dksi(ksi,eta,zeta),funDN8Dksi(ksi,eta,zeta);
                funDN1Deta(ksi,eta,zeta),funDN2Deta(ksi,eta,zeta),funDN3Deta(ksi,eta,zeta),funDN4Deta(ksi,eta,zeta), ...
                funDN5Deta(ksi,eta,zeta),funDN6Deta(ksi,eta,zeta),funDN7Deta(ksi,eta,zeta),funDN8Deta(ksi,eta,zeta);
                funDN1Dzeta(ksi,eta,zeta),funDN2Dzeta(ksi,eta,zeta),funDN3Dzeta(ksi,eta,zeta),funDN4Dzeta(ksi,eta,zeta), ...
                funDN5Dzeta(ksi,eta,zeta),funDN6Dzeta(ksi,eta,zeta),funDN7Dzeta(ksi,eta,zeta),funDN8Dzeta(ksi,eta,zeta)] * ...
                [pt1x,pt1y,pt1z;pt2x,pt2y,pt2z;pt3x,pt3y,pt3z;pt4x,pt4y,pt4z;pt5x,pt5y,pt5z;pt6x,pt6y,pt6z;pt7x,pt7y,pt7z;pt8x,pt8y,pt8z];
            
            Jacobian = det(J);
            InvJ = inv(J);
            
            % ------ Compute [DN] matrix ------
            DNOrig = [InvJ zeros(3,3) zeros(3,3); zeros(3,3) InvJ zeros(3,3); zeros(3,3) zeros(3,3) InvJ] * ...
                [funDN1Dksi(ksi,eta,zeta) 0 0 funDN2Dksi(ksi,eta,zeta) 0 0 funDN3Dksi(ksi,eta,zeta) 0 0 funDN4Dksi(ksi,eta,zeta) 0 0 ...
                funDN5Dksi(ksi,eta,zeta) 0 0 funDN6Dksi(ksi,eta,zeta) 0 0 funDN7Dksi(ksi,eta,zeta) 0 0 funDN8Dksi(ksi,eta,zeta) 0 0;
                funDN1Deta(ksi,eta,zeta) 0 0 funDN2Deta(ksi,eta,zeta) 0 0 funDN3Deta(ksi,eta,zeta) 0 0 funDN4Deta(ksi,eta,zeta) 0 0 ...
                funDN5Deta(ksi,eta,zeta) 0 0 funDN6Deta(ksi,eta,zeta) 0 0 funDN7Deta(ksi,eta,zeta) 0 0 funDN8Deta(ksi,eta,zeta) 0 0;
                funDN1Dzeta(ksi,eta,zeta) 0 0 funDN2Dzeta(ksi,eta,zeta) 0 0 funDN3Dzeta(ksi,eta,zeta) 0 0 funDN4Dzeta(ksi,eta,zeta) 0 0 ...
                funDN5Dzeta(ksi,eta,zeta) 0 0 funDN6Dzeta(ksi,eta,zeta) 0 0 funDN7Dzeta(ksi,eta,zeta) 0 0 funDN8Dzeta(ksi,eta,zeta) 0 0;
                0 funDN1Dksi(ksi,eta,zeta) 0 0 funDN2Dksi(ksi,eta,zeta) 0 0 funDN3Dksi(ksi,eta,zeta) 0 0 funDN4Dksi(ksi,eta,zeta) 0 ...
                0 funDN5Dksi(ksi,eta,zeta) 0 0 funDN6Dksi(ksi,eta,zeta) 0 0 funDN7Dksi(ksi,eta,zeta) 0 0 funDN8Dksi(ksi,eta,zeta) 0 ;
                0 funDN1Deta(ksi,eta,zeta) 0 0 funDN2Deta(ksi,eta,zeta) 0 0 funDN3Deta(ksi,eta,zeta) 0 0 funDN4Deta(ksi,eta,zeta) 0 ...
                0 funDN5Deta(ksi,eta,zeta) 0 0 funDN6Deta(ksi,eta,zeta) 0 0 funDN7Deta(ksi,eta,zeta) 0 0 funDN8Deta(ksi,eta,zeta) 0 ;
                0 funDN1Dzeta(ksi,eta,zeta) 0 0 funDN2Dzeta(ksi,eta,zeta) 0 0 funDN3Dzeta(ksi,eta,zeta) 0 0 funDN4Dzeta(ksi,eta,zeta) 0 ...
                0 funDN5Dzeta(ksi,eta,zeta) 0 0 funDN6Dzeta(ksi,eta,zeta) 0 0 funDN7Dzeta(ksi,eta,zeta) 0 0 funDN8Dzeta(ksi,eta,zeta) 0 ;
                0 0 funDN1Dksi(ksi,eta,zeta) 0 0 funDN2Dksi(ksi,eta,zeta) 0 0 funDN3Dksi(ksi,eta,zeta) 0 0 funDN4Dksi(ksi,eta,zeta) ...
                0 0 funDN5Dksi(ksi,eta,zeta) 0 0 funDN6Dksi(ksi,eta,zeta) 0 0 funDN7Dksi(ksi,eta,zeta) 0 0 funDN8Dksi(ksi,eta,zeta) ;
                0 0 funDN1Deta(ksi,eta,zeta) 0 0 funDN2Deta(ksi,eta,zeta) 0 0 funDN3Deta(ksi,eta,zeta) 0 0 funDN4Deta(ksi,eta,zeta) ...
                0 0 funDN5Deta(ksi,eta,zeta) 0 0 funDN6Deta(ksi,eta,zeta) 0 0 funDN7Deta(ksi,eta,zeta) 0 0 funDN8Deta(ksi,eta,zeta) ;
                0 0 funDN1Dzeta(ksi,eta,zeta) 0 0 funDN2Dzeta(ksi,eta,zeta) 0 0 funDN3Dzeta(ksi,eta,zeta) 0 0 funDN4Dzeta(ksi,eta,zeta) ...
                0 0 funDN5Dzeta(ksi,eta,zeta) 0 0 funDN6Dzeta(ksi,eta,zeta) 0 0 funDN7Dzeta(ksi,eta,zeta) 0 0 funDN8Dzeta(ksi,eta,zeta)];
            DN = DNOrig;
            
            % ------ Here approximate Dg(x+u)=Df(x) ------
            DfEle = [DfDx(ptOfx-DfDxStartx, ptOfy-DfDxStarty, ptOfz-DfDxStartz);
                DfDy(ptOfx-DfDxStartx, ptOfy-DfDxStarty, ptOfz-DfDxStartz);
                DfDz(ptOfx-DfDxStartx, ptOfy-DfDxStarty, ptOfz-DfDxStartz)];
            
            % ------ Only assemble stiffness in the first step ------
            if (stepwithinwhile==1)
                %A(tempIndexU,tempIndexU) = A(tempIndexU,tempIndexU) + (N'*DfEle)*(N'*DfEle)' + alpha*(DN')*DN;
                tempA = tempA + (N'*DfEle)*(N'*DfEle)' + alpha*(DN')*DN;
            end
            
            % ------ Construct b vector ------
            %temp1 = [ptOfx;ptOfy;ptOfz] + N*U(tempIndexU);
            %temp2 = ((Img1(ptOfx,ptOfy,ptOfz) - fungInterpolation_g3(temp1(1),temp1(2),temp1(3),Img2(floor(temp1(1))-1:floor(temp1(1))+2, ...
            %        floor(temp1(2))-1:floor(temp1(2))+2, floor(temp1(3))-1:floor(temp1(3))+2))) * (N'*DfEle));
            %b(tempIndexU) = b(tempIndexU) + tempb - (alpha*(DN')*DN)*U(tempIndexU);
            temp2 = ((Img1(ptOfx,ptOfy,ptOfz) - tempg(tempjj)) * (N'*DfEle));
            tempb = tempb + temp2 - (alpha*(DN')*DN)*U(tempIndexU);
            
            %end   % for ptOfx = pt1x:pt7x
            %end    % for ptOfy = pt1y:pt7y
            %end   % for ptOfz = pt1z:pt7z
        end
        
        
        % --- Store tempA and tempb ---
        if (stepwithinwhile==1)
            [IndexAXX,IndexAYY] = ndgrid(tempIndexU,tempIndexU);
            INDEXAI = [INDEXAI;IndexAXX(:)]; INDEXAJ = [INDEXAJ;IndexAYY(:)]; INDEXAVAL = [INDEXAVAL;tempA(:)]; %INDEXAREG = [INDEXAREG;tempAreg(:)];
        end
        INDEXBI = [INDEXBI;tempIndexU(:)]; INDEXBVAL = [INDEXBVAL;tempb(:)];
        
    end
    
    close(hbar);
    
    if (stepwithinwhile==1)
        A = sparse(INDEXAI,INDEXAJ,INDEXAVAL,FEMSize,FEMSize);
    end
    b = sparse(INDEXBI,ones(length(INDEXBI),1),INDEXBVAL,FEMSize,1);
    
    %%%%%%%%% Please ignore these lines for considering boundary conditions %%%%%%%%
    % coordsIndexInvolved = unique([0;elementsFEM(:)]); % Need modification for triangle elementsFEM
    %
    % % In adaptive mesh, using the following code:
    % UIndexInvolved = zeros(DIM*(length(coordsIndexInvolved)-1),1);
    %
    % % Not including the first 0-th entry
    % for tempi = 1:(size(coordsIndexInvolved,1)-1)
    %     UIndexInvolved(3*tempi-2:3*tempi) = [3*coordsIndexInvolved(tempi+1)-2; ...
    %                 3*coordsIndexInvolved(tempi+1)-1; 3*coordsIndexInvolved(tempi+1)];
    % end
    %
    % % ========= Set Dirichlet and Neumann boundary conditions =========
    % if isempty(dirichlet) ~= 1
    %     dirichlettemp = [3*dirichlet(:); 3*dirichlet(:)-1; 3*dirichlet(:)-2];
    % else
    %     dirichlettemp = [];
    % end
    % if isempty(neumann) ~= 1
    %     neumanntemp = [3*neumann(:,1); 3*neumann(:,1)-1; 3*neumann(:,1)-2; 3*neumann(:,2); 3*neumann(:,2)-1; 3*neumann(:,2)-2];
    % else
    %     neumanntemp = [];
    % end
    % FreeNodes = setdiff(UIndexInvolved,unique([dirichlettemp]));
    
    % ========= Neumann conditions ===========
    % Last step boundary condition force
    % BCForce = -Global_NodalStrainAvg(coordinatesFEM,elementsFEM,Uhat);
    % for tempj = 1:size(neumann,1)
    %     b(2*neumann(tempj,1:2)-1) = b(2*neumann(tempj,1:2)-1) + norm(coordinatesFEM(neumann(tempj,1),:)-coordinatesFEM(neumann(tempj,2),:)) ...
    %         * ( BCForce(4*neumann(tempj,1:2)-3) * neumann(tempj,3) + BCForce(4*neumann(tempj,1:2)-1) * neumann(tempj,4) );
    %     b(2*neumann(tempj,1:2))   = b(2*neumann(tempj,1:2)) + norm(coordinatesFEM(neumann(tempj,1),:)-coordinatesFEM(neumann(tempj,2),:)) ...
    %         * ( BCForce(4*neumann(tempj,1:2)-1) * neumann(tempj,3) + BCForce(4*neumann(tempj,1:2)) * neumann(tempj,4) );
    % end
    
    % % ========= Dirichlet conditions ==========
    % UhatOld = Uhat;
    % UhatNew = sparse(DIM*FEMSize + NodesPerEle*DIM, 1);
    % UhatNew(3*unique(dirichlet)) = U(3*unique(dirichlet));
    % UhatNew(3*unique(dirichlet)-1) = U(3*unique(dirichlet)-1);
    % UhatNew(3*unique(dirichlet)-2) = U(3*unique(dirichlet)-2);
    
    % ========= Solve FEM problem ===========
    W = A\b;
    
    normW = norm(W)/sqrt(size(W,1));
    normOfW(stepwithinwhile) = normW;
    timeICGN(stepwithinwhile) = toc; 
    U = reshape(U,length(U),1); W = reshape(W,length(W),1);
    
    disp(['normW = ',num2str(normW),' at iter ',num2str(stepwithinwhile),'; time cost = ',num2str(toc),'s']);
    
    if stepwithinwhile == 1
        normWOld = normW*10;
    else
        normWOld = normOfW(stepwithinwhile-1);
    end
    
    if (normW < tol) || ((normW/normWOld > 0.9) && (normW/normWOld < 1))
        U = U + W;
        break;
    elseif (normW >= tol  && normW < (0.1/tol)) % || ((normW/normWOld >= 1) && (normW/normWOld < 100)))
        U = U + W;
    else
        warning('Get diverged in Global_ICGN!!!')
        break;
    end
    
    
    
end

TotalTimeICGN = sum(timeICGN);
disp(['Elapsed time is ',num2str(TotalTimeICGN),' seconds.']);

end



%% ========= subroutines for  FEM shape function derivatives ========
function a = funDN1Dksi(ksi,eta,zeta)
a = 1/8*(-1)*(1-eta)*(1-zeta);
end
function a = funDN2Dksi(ksi,eta,zeta)
a = 1/8*( 1)*(1-eta)*(1-zeta);
end
function a = funDN3Dksi(ksi,eta,zeta)
a = 1/8*( 1)*(1+eta)*(1-zeta);
end
function a = funDN4Dksi(ksi,eta,zeta)
a = 1/8*(-1)*(1+eta)*(1-zeta);
end
function a = funDN5Dksi(ksi,eta,zeta)
a = 1/8*(-1)*(1-eta)*(1+zeta);
end
function a = funDN6Dksi(ksi,eta,zeta)
a = 1/8*( 1)*(1-eta)*(1+zeta);
end
function a = funDN7Dksi(ksi,eta,zeta)
a = 1/8*( 1)*(1+eta)*(1+zeta);
end
function a = funDN8Dksi(ksi,eta,zeta)
a = 1/8*(-1)*(1+eta)*(1+zeta);
end

% ----------------------------------------------------
function a = funDN1Deta(ksi,eta,zeta)
a = 1/8*(1-ksi)*(-1)*(1-zeta);
end
function a = funDN2Deta(ksi,eta,zeta)
a = 1/8*(1+ksi)*(-1)*(1-zeta);
end
function a = funDN3Deta(ksi,eta,zeta)
a = 1/8*(1+ksi)*( 1)*(1-zeta);
end
function a = funDN4Deta(ksi,eta,zeta)
a = 1/8*(1-ksi)*( 1)*(1-zeta);
end
function a = funDN5Deta(ksi,eta,zeta)
a = 1/8*(1-ksi)*(-1)*(1+zeta);
end
function a = funDN6Deta(ksi,eta,zeta)
a = 1/8*(1+ksi)*(-1)*(1+zeta);
end
function a = funDN7Deta(ksi,eta,zeta)
a = 1/8*(1+ksi)*( 1)*(1+zeta);
end
function a = funDN8Deta(ksi,eta,zeta)
a = 1/8*(1-ksi)*( 1)*(1+zeta);
end

% ----------------------------------------------------
function a = funDN1Dzeta(ksi,eta,zeta)
a = 1/8*(1-ksi)*(1-eta)*(-1);
end
function a = funDN2Dzeta(ksi,eta,zeta)
a = 1/8*(1+ksi)*(1-eta)*(-1);
end
function a = funDN3Dzeta(ksi,eta,zeta)
a = 1/8*(1+ksi)*(1+eta)*(-1);
end
function a = funDN4Dzeta(ksi,eta,zeta)
a = 1/8*(1-ksi)*(1+eta)*(-1);
end
function a = funDN5Dzeta(ksi,eta,zeta)
a = 1/8*(1-ksi)*(1-eta)*( 1);
end
function a = funDN6Dzeta(ksi,eta,zeta)
a = 1/8*(1+ksi)*(1-eta)*( 1);
end
function a = funDN7Dzeta(ksi,eta,zeta)
a = 1/8*(1+ksi)*(1+eta)*( 1);
end
function a = funDN8Dzeta(ksi,eta,zeta)
a = 1/8*(1-ksi)*(1+eta)*( 1);
end


