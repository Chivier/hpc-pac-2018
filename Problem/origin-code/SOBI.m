function W = SOBI(EEG_DATA, TAU, fs, jthresh)

    % The following code is prepared by The Laboratory of Neuroscience for
    % Education (NfE lab) at The University of Hong Kong for the use in PAC2018
    % competition only.
    %
    % Please note that due to the time limit, the code has not been cleaned up,
    % optimized and structrued. Therefore, the comprehensibility of the code
    % can be low.
    %
    % NfE lab owns the copyright of this program. The program is not allowed to be redistributed.

    if TAU(1) ~= 0; TAU = [0, TAU];  end

    for section = 1:1

        lblk = 1000;

        ms_tau = TAU; hertz = fs;

        if nargin == 0
            error('Warning: Must have input argument!');
        elseif nargin == 1
            msg = ('An A/D Hertz value was not entered, 1000Hz was assumed');
            disp(msg);
            hertz = 1000;
        end

        [discard, ms_tau_size] = size(ms_tau);

        for check1 = 1:ms_tau_size
            ms_tau(check1)

            if ms_tau(check1) < ceil(1000 / hertz)

                if ms_tau(check1) == 0
                else
                    msg = sprintf('Millisecond value for tau is too small use values larger than %d, the value detected is %d', ceil(1000 / hertz), ms_tau(check1));
                    disp(msg);
                end

            end

        end

        sec_tau = ms_tau / 1000;
        sample_tau = sec_tau * hertz;
        rounded_sample_tau = round(sample_tau);

        for check2 = 1:ms_tau_size

            if sample_tau < 1
                msg = sprintf('A tau value of 0 has been detected for millisecond value %d, increase value and rerun program if not desired.', ms_tau(check2));
                disp(msg);
            end

        end

        for check3 = 1:ms_tau_size

            for repeat = check3:ms_tau_size

                if check3 == repeat
                elseif rounded_sample_tau(check3) == -1
                elseif rounded_sample_tau(check3) == rounded_sample_tau(repeat)
                    msg = sprintf('A repeated value has been detected and deleted from the tau array, this value is %d', rounded_sample_tau(check3));
                    disp(msg);
                    rounded_sample_tau(repeat) = -1;
                end

            end

        end

        array_count = 0;
        opt_tau = 0;

        for make_tau = 1:ms_tau_size

            if rounded_sample_tau(make_tau) == -1
            else
                array_count = array_count + 1;
                opt_tau(array_count) = rounded_sample_tau(make_tau);
            end

        end

        TAU = opt_tau;

        Ntau = length(TAU);
        taumax = max(TAU);

        tcpu0 = cputime;
        Nsn = size(EEG_DATA, 1);
        Rtau_all = zeros(Nsn, Nsn, Ntau);

        t_end = 1 + size(EEG_DATA, 2) - lblk - taumax;

        TT = 0;

        for t0 = 1:lblk:t_end

            [XEEG] = EEG_DATA(:, t0:t0 + lblk + taumax - 1);

            X = XEEG - mean(XEEG, 2) * ones(1, lblk + taumax);

            for i = 1:Ntau,
                taui = TAU(i);
                RR = (X(:, 1:lblk) * X(:, 1 + taui:lblk + taui)');
                Rtau_all(:, :, i) = Rtau_all(:, :, i) + 0.5 * (RR + RR');
            end

            TT = TT + lblk;

            fprintf('Current block samples: %d-%d  Cpu time: %g\n', ...
                t0 + 1, t0 + lblk, cputime - tcpu0);

        end

        Rtau_all = (1 / TT) * Rtau_all;

    end

    for section = 2:2

        cputime0 = cputime;

        [Nsn, Nsn, Ntau] = size(Rtau_all); Ntau = Ntau - 1;
        N_princComp = Nsn;

        if nargin < 3, N_princComp = Nsn;  end
        if nargin < 4, method = 0;  end

        JOINT_DIAG = 0;
        TREE = 0;
        JOINT_DIAG = 1;

        R0 = Rtau_all(:, :, 1);

        [Nsn, ttt] = size(R0);

        N = Nsn;

        [S0, lam0] = eig(0.5 * (R0 + R0')); lam0 = diag(lam0);
        [lm0, I0] = sort(lam0);
        II0 = I0(Nsn:-1:Nsn - N + 1);
        lam1 = lam0(II0);
        S1 = S0(:, II0);
        B = (S1 * diag(sqrt(1 ./ (lam1 + 1e-20))))';
        lambda = lm0(Nsn:-1:1);

        N = N_princComp;
        B = B(1:N, :);

        disp('Presphering correlation matrices..');
        Rtau_presphered = zeros(N, N, Ntau);

        for i = 1:Ntau
            Rtau_presphered(:, :, i) = B * Rtau_all(:, :, i + 1) * B';
        end

        if JOINT_DIAG
            RRR = zeros(N, N * Ntau);
            RRR(:) = Rtau_presphered(:);

            A = RRR; B_mult = B;

            [m, nm] = size(A);

            B = [1 ,0, 0; 0 ,1, 1; 0 ,- 1, 1];
            Bt = B';
            Ip = zeros(1, nm);
            Iq = zeros(1, nm);
            g = zeros(3, nm);
            g = zeros(3, m);
            G = zeros(2, 2);
            vcp = zeros(3, 3);
            D = zeros(3, 3);
            la = zeros(3, 1);
            K = zeros(3, 3);
            angles = zeros(3, 1);
            pair = zeros(1, 2);
            G = zeros(3);
            c = 0;
            s = 0;

            for j = 1:8 eval(['flag' num2str(j) '=1;']);  end

            V = eye(m);
            encore = 1;
            iter_count = 0;

            while encore,
                encore = 0;
                smax = 0;

                for p = 1:m - 1, Ip = p:m:nm;

                    for q = p + 1:m, Iq = q:m:nm;

                        g = [A(p, Ip) - A(q, Iq); A(p, Iq); A(q, Ip)];
                        [vcp, D] = eig(real(B * (g * g') * Bt));
                        [la, K] = sort(diag(D));
                        angles = vcp(:, K(3));

                        if angles(1) < 0, angles = -angles; end;
                            c = sqrt(0.5 + angles(1) / 2);
                            s = 0.5 * (angles(2) - j * angles(3)) / c;

                            smax = max(s, smax);

                            if abs(s) > jthresh,
                                encore = 1;
                                pair = [p; q];
                                G = [c, - conj(s); s c];
                                V(:, pair) = V(:, pair) * G;
                                A(pair, :) = G' * A(pair, :);
                                A(:, [Ip Iq]) = [c * A(:, Ip) + s * A(:, Iq) ,- conj(s) * A(:, Ip) + c * A(:, Iq)];

                            end

                        end

                    end

                    W_unscaled = V' * B_mult;

                    unscaled_W = W_unscaled;

                    [m, n] = size(unscaled_W);
                    W_scaled = zeros(m, n);

                    for j = 1:m
                        row = unscaled_W(j, :).^2;
                        scaling_factor = sqrt(sum(row));
                        W_scaled(j, :) = unscaled_W(j, :) ./ scaling_factor;
                    end

                    W_scaled = real(W_scaled);

                    iter_count = iter_count + 1;
                    fprintf('Iteration %d accuracy: %e\n', iter_count, smax);

                end

                D = A;

                Wjoint = V; Djoint = D;
                W_xu = Wjoint' * B_mult;

                cputime3 = cputime - cputime0;

            end

            unscaled_W = W_xu;
            [m, n] = size(unscaled_W);
            W_scaled = zeros(m, n);

            for i = 1:m
                row = unscaled_W(i, :).^2;
                scaling_factor = sqrt(sum(row));
                W_scaled(i, :) = unscaled_W(i, :) ./ scaling_factor;
            end

            W_scaled = real(W_scaled);

            W = W_scaled;

        end
