% ============================================================
% Laplace equation solver (2D) by finite difference + Gauss-Seidel
% Domain: [0,1]x[0,1]
% BC:
%   Left  : phi = 1  (Dirichlet)
%   Right : phi = 0  (Dirichlet)
%   Top/Bottom : dphi/dy = 0 (Neumann)
% ============================================================
clear; clc;close all;


% grid
Nx = 51; Ny = 51;  %xとyの分割数
Lx = 1.0; Ly = 1.0; %ドメインの大きさ
dx = Lx/(Nx-1); dy = Ly/(Ny-1); %格子間隔
x = linspace(0, Lx, Nx);  y = linspace(0, Ly, Ny); %x,y座標
phi = zeros(Nx, Ny); %ポテンシャルの初期化

% BCの設定
%BCの設定とは、境界条件を指定することです。ここでは、左端と右端のDirichlet境界条件を設定しています。左端ではポテンシャルphiが1に固定され、右端ではポテンシャルphiが0に固定されます。これにより、流れの方向が左から右になることが示されます。上端と下端はNeumann境界条件で、dphi/dy = 0となっているため、これらの境界ではポテンシャルの変化がないことを意味します。
%phi(1, :) = 1; %左端のDirichlet BC  
%phi(Nx, :) = 0; %右端のDirichlet BC

% 上下境界条件をNeumannではなく，Dirichletにする場合は、以下のように設定します。
%phi(:, 1) = 0; %下端のDirichlet BC
%phi(:, Ny) = 0; %上端のDirichlet BC

%境界条件をランダム関数を用いても同じ解が得られることを確認するために、以下のように設定します。
%phi = rand(Nx, Ny); %左端のDirich

% ソルバーを設定する
max_iter = 10000; tol = 1e-8; %最大反復回数と収束判定の閾値

% iterative solver (Gauss-Seidel)繰り返し計算という意味
for iter = 1:max_iter
    phi_old = phi; %前の反復のポテンシャルを保存

        %左右のDirichlet BCの適用(BCとはBoundary Conditionの略で、境界条件のことを指します。ここでは、左右のDirichlet BCを適用しています。Dirichlet BCは、特定の値にポテンシャルを固定する境界条件です。左端ではphiを1に固定し、右端ではphiを0に固定しています。これにより、流れの方向が左から右になることが示されます。)
    phi(1, :) = 1; %左端のDirichlet BC．何をしているのかというと，左端のDirichlet BCを適用しています。phi(1, :) = 1というコードは、左端の全ての点に対してポテンシャルphiを1に固定することを意味します。これにより、流れの方向が左から右になることが示されます。同様に、右端のDirichlet BCも適用され、phi(Nx, :) = 0となります。
    phi(Nx, :) = 0; %右端のDirichlet BC．何をしているのかというと，右端のDirichlet BCを適用しています。phi(Nx, :) = 0というコードは、右端の全ての点に対してポテンシャルphiを0に固定することを意味します。これにより、流れの方向が左から右になることが示されます。同様に、左端のDirichlet BCも適用され、phi(1, :) = 1となります。

    % Neumann BCの適用（上端と下端）
    phi(:, 1) = phi(:, 2);   %何をしているのかというと，下端のNeumann BCを適用しています。dphi/dy = 0という条件は、ポテンシャルがy方向に変化しないことを意味します。したがって、下端のポテンシャルphi(:, 1)は、そのすぐ上の点phi(:, 2)と同じ値になります。同様に、上端のNeumann BCも同様に適用され、phi(:, Ny)はそのすぐ下の点phi(:, Ny-1)と同じ値になります。
    phi(:, Ny) = phi(:, Ny-1); %何をしているのかというと，上端のNeumann BCを適用しています。dphi/dy = 0という条件は、ポテンシャルがy方向に変化しないことを意味します。したがって、上端のポテンシャルphi(:, Ny)は、そのすぐ下の点phi(:, Ny-1)と同じ値になります。同様に、下端のNeumann BCも同様に適用され、phi(:, 1)はそのすぐ上の点phi(:, 2)と同じ値になります。
    for i = 2:Nx-1
        for j = 2:Ny-1
            % ラプラス方程式の離散化
            phi(i, j) = 0.25 * (phi(i+1, j) + phi(i-1, j) + phi(i, j+1) + phi(i, j-1));% 4点平均法
        end
    end
    


    
    % 収束判定
    if max(max(abs(phi - phi_old))) < tol
        fprintf('Converged after %d iterations.\n', iter);
        break;
    end
end

% 結果のプロット
figure;
contourf(x, y, phi', 40); %等高線プロット
colorbar; xlabel('x'); ylabel('y'); title('Potential Flow (Laplace Equation)');
%色が何を示しているのかというと、等高線プロットにおいて、色はポテンシャルphiの値を表しています。色のグラデーションは、ポテンシャルの大きさに対応しており、通常は高い値が暖色系（赤や黄色）で、低い値が寒色系（青や緑）で表示されます。これにより、流れのパターンやポテンシャルの分布を視覚的に理解することができます。
%色が暖色であれば，ポテンシャルが高いことを示し、寒色であれば、ポテンシャルが低いことを示しています。つまり，今回であれば流れは左から右に向かっていることがわかります。

fprintf('Converged after %d iterations.\n', iter);