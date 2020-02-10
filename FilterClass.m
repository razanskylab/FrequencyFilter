classdef FilterClass < BaseClass

  properties
    df(1, 1) {mustBeNumeric, mustBeFinite} = 250; % sampling rate

    freq(1, :) {mustBeNumeric, mustBeFinite} = [2 40];
    % freq = [1e6];
    % [] = no filtering
    % [x] = high pass
    % [x y] = band pass

    order(1, :) {mustBeNumeric} = 1;
    % higher orders work as well, are slower but result in better
    % signal representation

    filtType(1, 1) {mustBeNumeric, mustBeFinite} = 1; %
    % 1 - butterworth (butter)
    % 2 - Chebyshev Type I (cheby1)

    tech(1, 1) {mustBeNumeric, mustBeFinite} = 1;
    % 0 - use matlab filtfilt
    % 1 - use file exchange filtfiltM - faster and accepts single

    filtMode = '2d'; % '2d' | '3d'
  end

  properties (SetAccess = private)
    % discrete-time zero-pole-gain representation of a given digital filter
    % z;
    % p;
    % k;
    % % second-order sections form
    % sos;
    % g;
    % Transfer function coefficients, specified as vectors.
    % NOTE Should not be used, as z,p,k are more stable!
    b;
    a;
  end

  properties (Dependent = true)
    normF; % normalized frequency, the one used in most filter definitions
    % [sos,g] = zp2sos(z,p,k);
  end

  properties
    % defined in base class but this way we can have set/get in sublcasses, which
    % is % needed for FOAM processor
    silent(1, 1) {mustBeNumericOrLogical} = false;
    verboseOutput(1, 1) {mustBeNumericOrLogical} = false; % more detailed output to workspace...
    verbosePlotting(1, 1) {mustBeNumericOrLogical} = false; % more figures...
    figureVisibility(1, :) char {mustBeMember(figureVisibility, {'on', 'off'})} = 'on';
  end

  properties (Constant)
    PASS_BAND_RIPPLE = 0.25; % [dB], less means flatter filter response
    % Peak-to-peak passband ripple, in decibels
  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  methods
    % define filters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function Define(FC)
      tic;
      FC.VPrintF('   Defining filter specs, using:\n')
      % get strings expressing filter bandwidth or cutoff freq. for high pass
      % filter
      if numel(FC.freq) == 2
        FC.VPrintF('      Bandpass (%3.1f MHz to %3.1f MHz)\n', ...
          FC.freq(1), FC.freq(2));
      else
        FC.VPrintF('      Highpass (passband = %3.1f MHz)\n', FC.freq(1));
      end

      switch FC.filtType
          % butterworth filter -----------------------------------------------------
        case 1
          % 1 - first-order butterworth, not good for strong dc bias
          if length(FC.freq) == 1
            [FC.b, FC.a] = butter(FC.order, FC.normF, 'high');
          elseif length(FC.freq) == 2
            [FC.b, FC.a] = butter(FC.order, FC.normF, 'bandpass');
          end

          FC.VPrintF('      Butterworth filter (order %i)\n', FC.order);
          % Chebyshev Type I filter ------------------------------------------------
        case 2
          % 2 - first-order Chebyshev Type I, sharper and better for DC bias
          passBandRipple = 5; %Peak-to-peak passband ripple in db

          if length(FC.freq) == 1
            [FC.b, FC.a] = cheby1(FC.order, passBandRipple, FC.normF, 'high');
          elseif length(FC.freq) == 2
            [FC.b, FC.a] = cheby1(FC.order, passBandRipple, FC.normF, 'bandpass');
          end

          FC.VPrintF('      Chebyshev Type I filter (order %i)\n', FC.order);
        otherwise
          error('Unknown filter type, FC.filtType must be 1 or 2!')
      end

      % Convert zero-pole-gain filter parameters to second-order sections form
      % [FC.sos,FC.g] = zp2sos(FC.z,FC.p,FC.k);
      % Convert zero-pole-gain filter parameters to transfer function form
      % [FC.b,FC.a] = zp2tf(FC.z,FC.p,FC.k);
      FC.VPrintF('      Done in %2.2f s.\n', toc);
    end

    % apply filters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [usShots] = Apply(FC, usShots)

      if nargin == 1
        error('Filter.Apply requires raw-data as input argument!')
      end

      if any(isempty([FC.b FC.a]))
        FC.VPrintF('   Defining filter before applying...');
        FC.Define();
      end

      tic();
      dataType = class(usShots);

      if FC.tech == 0% use matlab filtfilt which requires double and is slower
        % FC.VPrintF('   Filtering using filtfilt...')
        usShots = filtfilt(FC.b, FC.a, double(usShots));
      end

      if FC.tech == 1% using faster FiltM funciton form Maltab file exchange
        % FC.VPrintF('   Filtering using FiltFiltM...')
        usShots = FiltFiltM(FC.b, FC.a, single(usShots));
      end

      usShots = cast(usShots, dataType); % restore data type to what it was before
      FC.VPrintF('done in %2.2f s.\n', toc);
    end

    % apply filter to volume %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [volData] = Apply_Vol(FC, volData)
      % filters volumes of the form xyz!!! so filters last dimension!
      if nargin == 1
        error('Filter.Apply requires raw-data as input argument!')
      end

      if any(isempty([FC.b FC.a]))
        FC.VPrintF('   Defining filter before applying:\n');
        FC.Define();
      end

      tic();

      [nX, nY, nZ] = size(volData);
      dataType = class(volData);

      switch FC.filtMode
        case '2d'
          % reshape to 2d, then apply filter ---------------------------------------
          volData = reshape(volData, nX * nY, nZ)';
          % need to permute to filter along right dimension3

          if FC.tech == 0% use matlab filtfilt which requires double and is slower
            FC.VPrintF('   Filtering using filtfilt...')
            volData = filtfilt(FC.b, FC.a, double(volData));
          end

          if FC.tech == 1% using faster FiltM funciton form Maltab file exchange
            FC.VPrintF('   Filtering using FiltFiltM...')
            volData = FiltFiltM(FC.b, FC.a, volData);
          end

          volData = reshape(volData', nX, nY, nZ);
          volData = cast(volData, dataType);
        case '3d'
          % filters along first dimension, accepts single and double
          volData = filtfiltj(FC.b, FC.a, volData);
        otherwise
          short_warn('Unknow filter mode! Not doing anything!');
          volData = [];
      end

    end

    % Plot %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [] = Plot(FC, plotAx)
      % plot_filter_response(FC.df,FC.b,FC.a);
      if nargin == 1
        plotAx = gca;
      end

      % get info on filter response for later plotting
      nFiltPoints = 250; % the less the faster...
      [highPasAmpResp, ~] = freqz(FC.b, FC.a, nFiltPoints, FC.df);
      [highPasPhaseResp, filtFreq] = phasez(FC.b, FC.a, nFiltPoints, FC.df);
      highPasPhaseResp = rad2deg(highPasPhaseResp); % convert to deg, find that easier
      highPasAmpResp = highPasAmpResp * 100; % convert to%, also easier for me...

      yyaxis(plotAx, 'left')
      plot(plotAx, filtFreq, abs(highPasAmpResp));
      axis(plotAx, 'tight');
      ylabel(plotAx, 'Amplitude');
      ylim(plotAx, [0 110]);
      yyaxis(plotAx, 'right');
      plot(plotAx, filtFreq, abs(highPasPhaseResp));
      ylabel(plotAx, 'Phase Shift');
      axis(plotAx, 'tight');
      grid(plotAx, 'on');
      % xlim([0 180]);
      % title(sprintf('high pass filter (n = %i, f = %2.1f MHz)',FC.order,lowPassFreq*1e-6));
    end

    % Plot %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [] = Plot_Amp(FC, plotAx)
      % plot_filter_response(FC.df,FC.b,FC.a);
      if nargin == 1
        plotAx = gca;
      end

      % get info on filter response for later plotting
      nFiltPoints = 200; % the less the faster...
      [highPasAmpResp, filtFreq] = freqz(FC.b, FC.a, nFiltPoints, FC.df);
      highPasAmpResp = highPasAmpResp * 100; % convert to%, also easier for me...

      plot(plotAx, filtFreq, abs(highPasAmpResp));
      axis(plotAx, 'tight');
      ylabel(plotAx, 'Amplitude');
      axis(plotAx, 'tight');
      grid(plotAx, 'on');
      ylim(plotAx, [0 110]);
      % xlim([0 180]);
      % title(sprintf('high pass filter (n = %i, f = %2.1f MHz)',FC.order,lowPassFreq*1e-6));
    end

    % Debug_Stuf % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [] = Debug_Stuff(~, usShots)

      % define filter based on given filter parameters
      % then apply that filter using filtfilt or faster filtM
      if show_debug_plot()% plot here to not have to keep a copy of the raw us data
        figure();
        nPlots = 9; % plot this many signals
        [m, n] = find_subplot_dividers(nPlots);
        nShots = size(usShots, 1);
        plotSignalIdx = 1:round(nShots / nPlots):nShots;
        plotSignals = usShots(plotSignalIdx, :);

        for iPlot = 1:nPlots
          subplot(m, n, iPlot);
          plot(plotSignals(iPlot, :));
          axis tight;
        end

      end

      if show_debug_plot()% plot here to not have to keep a copy of the raw us data
        plotSignals = usShots(plotSignalIdx, :);

        for iPlot = 1:nPlots
          subplot(m, n, iPlot);
          hold on;
          plot(plotSignals(iPlot, :));
          axis tight;
        end

        sub_plot_title('Raw vs Filt Data');
      end

    end

    %%===========================================================================
    function normF = get.normF(FC)
      normF = FC.freq / (FC.df / 2);
    end

  end % end methods

end
