#include <TH1Plotter.h>

#include <TCanvas.h>
#include <TEfficiency.h>
#include <TF1.h>
#include <TFitResult.h>
#include <TLatex.h>
#include <TLine.h>
#include <TObject.h>
#include <TPave.h>
#include <TVirtualFitter.h>
#include <TGraphAsymmErrors.h>

#include <commandlinecfg.h>
#include <pool.h>
#include <utilities.h>

namespace plotIt {

    /*!
     * Compute the ratio between two histograms, taking into account asymetric error bars
     */
    std::shared_ptr<TGraphAsymmErrors> getRatio(TH1* a, TH1* b) {
        std::shared_ptr<TGraphAsymmErrors> g(new TGraphAsymmErrors(a));

        size_t npoint = 0;
        for (size_t i = 1; i <= (size_t) a->GetNbinsX(); i++) {
            float b1 = a->GetBinContent(i);
            float b2 = b->GetBinContent(i);

            if ((b1 == 0) || (b2 == 0))
                continue;

            float ratio = b1 / b2;

            float b1sq = b1 * b1;
            float b2sq = b2 * b2;

            float e1sq_up = a->GetBinErrorUp(i) * a->GetBinErrorUp(i);
            float e2sq_up = b->GetBinErrorUp(i) * b->GetBinErrorUp(i);

            float e1sq_low = a->GetBinErrorLow(i) * a->GetBinErrorLow(i);
            float e2sq_low = b->GetBinErrorLow(i) * b->GetBinErrorLow(i);

            float error_up = sqrt((e1sq_up * b2sq + e2sq_up * b1sq) / (b2sq * b2sq));
            float error_low = sqrt((e1sq_low * b2sq + e2sq_low * b1sq) / (b2sq * b2sq));

            //Set the point center and its errors
            g->SetPoint(npoint, a->GetBinCenter(i), ratio);
            //g->SetPointError(npoint, a->GetBinCenter(i) - a->GetBinLowEdge(i),
                    //a->GetBinLowEdge(i) - a->GetBinCenter(i) + a->GetBinWidth(i),
                    //error_low, error_up);
            g->SetPointError(npoint, 0, 0, error_low, error_up);
            npoint++;
        }

        g->Set(npoint);

        return g;
    }

  bool TH1Plotter::supports(TObject& object) {
    return object.InheritsFrom("TH1");
  }

  boost::optional<Summary> TH1Plotter::plot(TCanvas& c, Plot& plot) {
    c.cd();

    Summary global_summary;

    // Rescale and style histograms
    for (File& file: m_plotIt.getFiles()) {
      setHistogramStyle(file);

      TH1* h = dynamic_cast<TH1*>(file.object);

      if (file.type != DATA) {
        plot.is_rescaled = true;

        float factor = file.cross_section * file.branching_ratio / file.generated_events;

        if (! m_plotIt.getConfiguration().no_lumi_rescaling) {
          factor *= m_plotIt.getConfiguration().luminosity;
        }

        if (! CommandLineCfg::get().ignore_scales) {
          factor *= m_plotIt.getConfiguration().scale * file.scale;
        }

        h->Scale(factor);

        SummaryItem summary;
        summary.name = file.pretty_name;
        summary.process_id = file.id;

        double rescaled_integral_error = 0;
        double rescaled_integral = h->IntegralAndError(h->GetXaxis()->GetFirst(), h->GetXaxis()->GetLast(), rescaled_integral_error);

        summary.events = rescaled_integral;
        summary.events_uncertainty = rescaled_integral_error;

        // FIXME: Probably invalid in case of weights...

        // Bayesian efficiency
        // Taken from https://root.cern.ch/doc/master/TEfficiency_8cxx_source.html#l02428

        /*
        // Use a flat prior (equivalent to Beta(1, 1))
        float alpha = 1.;
        float beta = 1.;

        summary.efficiency = TEfficiency::BetaMean(integral + alpha, file.generated_events - integral + beta);
        summary.efficiency_uncertainty = TEfficiency::Bayesian(file.generated_events, integral, 0.68, alpha, beta, true) - summary.efficiency;
        */

        global_summary.add(file.type, summary);

        // Update all systematics for this file
        for (auto& syst: *file.systematics) {
          syst.update();

          syst.scale(factor);
          syst.rebin(plot.rebin);
        }

      } else {
        SummaryItem summary;
        summary.name = file.pretty_name;
        summary.process_id = file.id;
        summary.events = h->Integral();
        global_summary.add(file.type, summary);
      }

      h->Rebin(plot.rebin);

      // Add overflow to first and last bin if requested
      if (plot.show_overflow) {
        addOverflow(h, file.type, plot);

        if (file.type != DATA) {
            for (auto& syst: *file.systematics) {
                addOverflow(static_cast<TH1*>(syst.nominal_shape.get()), file.type, plot);
                addOverflow(static_cast<TH1*>(syst.up_shape.get()), file.type, plot);
                addOverflow(static_cast<TH1*>(syst.down_shape.get()), file.type, plot);
            }
        }
      }
    }

    // Build a THStack for MC files and a vector for signal
    float mcWeight = 0;
    std::shared_ptr<THStack> mc_stack;
    std::shared_ptr<TH1> mc_histo_stat_only;
    std::shared_ptr<TH1> mc_histo_syst_only;
    std::shared_ptr<TH1> mc_histo_stat_syst;

    std::shared_ptr<TH1> h_data;
    std::string data_drawing_options;

    std::vector<File> signal_files;

    // First pass. Create one histogram per group
    // Key is group name, value is group histogram
    std::vector<std::pair<std::string, std::shared_ptr<TH1>>> group_histograms;
    for (File& file: m_plotIt.getFiles()) {
        if (file.type != MC)
            continue;

        if (file.legend_group.empty())
            continue;

        auto it = std::find_if(group_histograms.begin(), group_histograms.end(), [&file](const std::pair<std::string, std::shared_ptr<TH1>>& item) {
            return item.first == file.legend_group;
        });

        TH1* nominal = dynamic_cast<TH1*>(file.object);

        // Do not bother with histograms with no entries
        if (nominal->GetEntries() == 0)
            continue;

        if (it == group_histograms.end()) {
            std::shared_ptr<TH1> h(dynamic_cast<TH1*>(nominal->Clone()));
            h->SetDirectory(nullptr);
            group_histograms.push_back(std::make_pair(file.legend_group, h));
        } else {
            it->second->Add(nominal);
        }
    }

    for (File& file: m_plotIt.getFiles()) {
      if (file.type == MC) {

        TH1* nominal = dynamic_cast<TH1*>(file.object);

        // Do not bother with histograms with no entries
        if (file.legend_group.empty() && nominal->GetEntries() == 0)
            continue;

        if (mc_stack.get() == nullptr)
          mc_stack = std::make_shared<THStack>("mc_stack", "mc_stack");

        auto it = std::find_if(group_histograms.begin(), group_histograms.end(), [&file](const std::pair<std::string, std::shared_ptr<TH1>>& item) {
            return item.first == file.legend_group;
        });

        if (!file.legend_group.empty() && it == group_histograms.end()) {
            // The group histogram has already been added to the stack, so
            // skip to the next one
            continue;
        } else if (it != group_histograms.end()) {
            auto n = it->second;
            TemporaryPool::get().add(n);

            nominal = n.get();

            // Since we will add this group histogram to the stack
            // remove it from the pool to avoid double addition
            group_histograms.erase(it);
        }

        mc_stack->Add(nominal, m_plotIt.getPlotStyle(file)->drawing_options.c_str());
        if (mc_histo_stat_only.get()) {
          mc_histo_stat_only->Add(nominal);
        } else {
          mc_histo_stat_only.reset( dynamic_cast<TH1*>(nominal->Clone()) );
          mc_histo_stat_only->SetDirectory(nullptr);
        }
        mcWeight += nominal->GetSumOfWeights();

      } else if (file.type == SIGNAL) {
        signal_files.push_back(file);
      } else if (file.type == DATA) {
        if (! h_data.get()) {
          h_data.reset(dynamic_cast<TH1*>(file.object->Clone()));
          h_data->SetDirectory(nullptr);
          h_data->Sumw2(false); // Disable SumW2 for data
          h_data->SetBinErrorOption((TH1::EBinErrorOpt) plot.errors_type);
          data_drawing_options += m_plotIt.getPlotStyle(file)->drawing_options;
        } else {
          h_data->Add(dynamic_cast<TH1*>(file.object));
        }
      }
    }

    if (plot.no_data || ((h_data.get()) && !h_data->GetSumOfWeights()))
      h_data.reset();

    if ((mc_histo_stat_only.get() && !mc_histo_stat_only->GetSumOfWeights())) {
      mc_histo_stat_only.reset();
      mc_stack.reset();
    }

    if (plot.normalized) {
      // Normalized each plot
      for (File& file: m_plotIt.getFiles()) {
        TH1* h = dynamic_cast<TH1*>(file.object);
        if (file.type == MC) {
          h->Scale(1. / fabs(mcWeight));
        } else if (file.type == SIGNAL) {
          h->Scale(1. / fabs(h->GetSumOfWeights()));
        }
      }

      if (h_data.get()) {
        h_data->Scale(1. / h_data->GetSumOfWeights());
      }
    }

    // Blind data if requested
    // It's not enough to put the bin content to zero, because
    // ROOT will show the marker, even with 'P'
    // The histogram is cloned, reset, and only the non-blinded bins are filled
    std::shared_ptr<TBox> m_blinded_area;
    if (!CommandLineCfg::get().unblind && h_data.get() && plot.blinded_range.valid()) {
        float start = plot.blinded_range.start;
        float end = plot.blinded_range.end;

        size_t start_bin = h_data->FindBin(start);
        size_t end_bin = h_data->FindBin(end);

        TH1* clone = static_cast<TH1*>(h_data->Clone());
        clone->SetDirectory(nullptr);

        h_data->Reset();
        h_data->Sumw2(false);

        for (size_t i = 0; i < start_bin; i++) {
            h_data->SetBinContent(i, clone->GetBinContent(i));
        }

        for (size_t i = end_bin + 1; i <= static_cast<size_t>(h_data->GetNbinsX()); i++) {
            h_data->SetBinContent(i, clone->GetBinContent(i));
        }

        delete clone;
    }

    if (mc_histo_stat_only.get()) {
      mc_histo_syst_only.reset(static_cast<TH1*>(mc_histo_stat_only->Clone()));
      mc_histo_syst_only->SetDirectory(nullptr);
      mc_histo_stat_syst.reset(static_cast<TH1*>(mc_histo_stat_only->Clone()));
      mc_histo_stat_syst->SetDirectory(nullptr);

      // Clear statistical errors
      for (uint32_t i = 1; i <= (uint32_t) mc_histo_syst_only->GetNbinsX(); i++) {
        mc_histo_syst_only->SetBinError(i, 0);
      }
    }

    if (mc_histo_syst_only.get() && plot.show_errors) {
      if (m_plotIt.getConfiguration().luminosity_error_percent > 0) {
        // Loop over all bins, and add lumi error
        for (uint32_t i = 1; i <= (uint32_t) mc_histo_syst_only->GetNbinsX(); i++) {
          float error = mc_histo_syst_only->GetBinError(i);
          float entries = mc_histo_syst_only->GetBinContent(i);
          float lumi_error = entries * m_plotIt.getConfiguration().luminosity_error_percent;

          mc_histo_syst_only->SetBinError(i, std::sqrt(error * error + lumi_error * lumi_error));
        }
      }

      // Key is systematics name, value is the combined systematics value for each bin
      std::map<std::string, std::vector<float>> combined_systematics_map;

      for (File& file: m_plotIt.getFiles()) {
        if (file.type == DATA || file.systematics->size() == 0)
          continue;

        for (auto& syst: *file.systematics) {

          std::vector<float>* combined_systematics;
          auto map_it = combined_systematics_map.find(syst.name());

          if (map_it == combined_systematics_map.end()) {
            combined_systematics = &combined_systematics_map[syst.name()];
            combined_systematics->resize(mc_histo_syst_only->GetNbinsX(), 0);
          } else {
            combined_systematics = &map_it->second;
          }

          TH1* nominal_shape = static_cast<TH1*>(syst.nominal_shape.get());
          TH1* up_shape = static_cast<TH1*>(syst.up_shape.get());
          TH1* down_shape = static_cast<TH1*>(syst.down_shape.get());

          if (! nominal_shape || ! up_shape || ! down_shape)
              continue;

          float total_syst_error = 0;
          // Systematics in each bin are fully correlated, as they come either from
          // a global variation, or for a shape variation. The total systematics error
          // is simply for sum of all errors in each bins
          //
          // However, we consider that different systematics in the same bin are totaly
          // uncorrelated. The total systematics errors is then the quadratic sum.
          for (uint32_t i = 1; i <= (uint32_t) mc_histo_syst_only->GetNbinsX(); i++) {
            float syst_error_up = std::abs(up_shape->GetBinContent(i) - nominal_shape->GetBinContent(i));
            float syst_error_down = std::abs(nominal_shape->GetBinContent(i) - down_shape->GetBinContent(i));

            // FIXME: Add support for asymetric errors
            float syst_error = std::max(syst_error_up, syst_error_down);

            total_syst_error += syst_error;

            // Only propagate uncertainties for MC, not signal
            if (file.type == MC)
                (*combined_systematics)[i - 1] += syst_error;
          }


          SummaryItem summary;
          summary.process_id = file.id;
          summary.name = syst.prettyName();
          summary.events_uncertainty = total_syst_error;

          global_summary.addSystematics(file.type, file.id, summary);
        }

      }

      // Combine all systematics in one
      // Consider that all the systematics are not correlated
      for (auto& combined_systematics: combined_systematics_map) {
        for (size_t i = 1; i <= (size_t) mc_histo_syst_only->GetNbinsX(); i++) {
          float total_error = mc_histo_syst_only->GetBinError(i);
          mc_histo_syst_only->SetBinError(i, std::sqrt(total_error * total_error + combined_systematics.second[i - 1] * combined_systematics.second[i - 1]));
        }
      }

      // Propagate syst errors to the stat + syst histogram
      for (uint32_t i = 1; i <= (uint32_t) mc_histo_syst_only->GetNbinsX(); i++) {
        float syst_error = mc_histo_syst_only->GetBinError(i);
        float stat_error = mc_histo_stat_only->GetBinError(i);
        mc_histo_stat_syst->SetBinError(i, std::sqrt(syst_error * syst_error + stat_error * stat_error));
      }
    }

    // Store all the histograms to draw, and find the one with the highest maximum
    std::vector<std::pair<TObject*, std::string>> toDraw = { std::make_pair(mc_stack.get(), ""), std::make_pair(h_data.get(), data_drawing_options) };
    for (File& signal: signal_files) {
      toDraw.push_back(std::make_pair(signal.object, m_plotIt.getPlotStyle(signal)->drawing_options));
    }

    // Remove NULL items
    toDraw.erase(
        std::remove_if(toDraw.begin(), toDraw.end(), [](std::pair<TObject*, std::string>& element) {
          return element.first == nullptr;
        }), toDraw.end()
      );

    if (!toDraw.size()) {
      std::cerr << "Error: nothing to draw." << std::endl;
      return boost::none;
    };

    // Sort object by minimum
    std::sort(toDraw.begin(), toDraw.end(), [&plot](std::pair<TObject*, std::string> a, std::pair<TObject*, std::string> b) {
        return (!plot.log_y) ? (getMinimum(a.first) < getMinimum(b.first)) : (getPositiveMinimum(a.first) < getPositiveMinimum(b.first));
      });

    float minimum = (!plot.log_y) ? getMinimum(toDraw[0].first) : getPositiveMinimum(toDraw[0].first);

    // Sort objects by maximum
    std::sort(toDraw.begin(), toDraw.end(), [](std::pair<TObject*, std::string> a, std::pair<TObject*, std::string> b) {
        return getMaximum(a.first) > getMaximum(b.first);
      });

    float maximum = getMaximum(toDraw[0].first);

    if ((!h_data.get() || !mc_histo_stat_only.get()))
      plot.show_ratio = false;

    std::shared_ptr<TPad> hi_pad;
    std::shared_ptr<TPad> low_pad;
    if (plot.show_ratio) {
      hi_pad = std::make_shared<TPad>("pad_hi", "", 0., 0.33333, 1, 1);
      hi_pad->Draw();
      hi_pad->SetTopMargin(TOP_MARGIN / .6666);
      hi_pad->SetLeftMargin(LEFT_MARGIN);
      hi_pad->SetBottomMargin(0.015);
      hi_pad->SetRightMargin(RIGHT_MARGIN);

      low_pad = std::make_shared<TPad>("pad_lo", "", 0., 0., 1, 0.33333);
      low_pad->Draw();
      low_pad->SetLeftMargin(LEFT_MARGIN);
      low_pad->SetTopMargin(1.);
      low_pad->SetBottomMargin(BOTTOM_MARGIN / .3333);
      low_pad->SetRightMargin(RIGHT_MARGIN);
      low_pad->SetTickx(1);

      hi_pad->cd();
      if (plot.log_y)
        hi_pad->SetLogy();

      if (plot.log_x) {
        hi_pad->SetLogx();
        low_pad->SetLogx();
      }
    }

    // Take into account systematics for maximum
    if (mc_histo_syst_only) {
      float maximum_with_errors = 0;
      for (size_t b = 1; b <= (size_t) mc_histo_stat_syst->GetNbinsX(); b++) {
        float local_maximum = mc_histo_stat_syst->GetBinContent(b) + mc_histo_stat_syst->GetBinErrorUp(b);
        maximum_with_errors = std::max(maximum_with_errors, local_maximum);
      }

      maximum = std::max(maximum, maximum_with_errors);
    }

    toDraw[0].first->Draw(toDraw[0].second.c_str());
    setRange(toDraw[0].first, plot.x_axis_range, plot.y_axis_range);

    float safe_margin = .20;
    if (plot.log_y)
      safe_margin = 8;

    if (! plot.y_axis_range.valid()) {
      maximum *= 1 + safe_margin;
      setMaximum(toDraw[0].first, maximum);

      if (minimum <= 0 && plot.log_y) {
        double old_minimum = minimum;
        minimum = 0.1;
        std::cout << "Warning: detected minimum is negative (" << old_minimum << ") but log scale is on. Setting minimum to " << minimum << std::endl;
      }

      if (!plot.log_y)
        minimum = minimum * (1 - std::copysign(safe_margin, minimum));

      if (plot.y_axis_show_zero && !plot.log_y)
        minimum = 0;

      setMinimum(toDraw[0].first, minimum);
    } else {
        maximum = plot.y_axis_range.end;
        minimum = plot.y_axis_range.start;
    }

    // First, draw MC
    if (mc_stack.get()) {
      mc_stack->Draw("same");

      // Clear all the possible stats box remaining
      mc_stack->GetHistogram()->SetStats(false);
      TIter next(mc_stack->GetHists());
      TH1* h = nullptr;
      while ((h = static_cast<TH1*>(next()))) {
          h->SetStats(false);
      }

      TemporaryPool::get().add(mc_stack);
    }

    // Then, if requested, errors
    if (mc_histo_stat_syst.get() && plot.show_errors) {
      mc_histo_stat_syst->SetMarkerSize(0);
      mc_histo_stat_syst->SetMarkerStyle(0);
      mc_histo_stat_syst->SetFillStyle(m_plotIt.getConfiguration().error_fill_style);
      mc_histo_stat_syst->SetFillColor(m_plotIt.getConfiguration().error_fill_color);

      mc_histo_stat_syst->Draw("E2 same");
      TemporaryPool::get().add(mc_histo_stat_syst);
    }

    // Then signal
    for (File& signal: signal_files) {
      std::string options = m_plotIt.getPlotStyle(signal)->drawing_options + " same";
      signal.object->Draw(options.c_str());
    }

    // And finally data
    if (h_data.get()) {
      data_drawing_options += " same";
      h_data->Draw(data_drawing_options.c_str());
      TemporaryPool::get().add(h_data);
    }

    // Set x and y axis titles, and default style
    for (auto& obj: toDraw) {
      setDefaultStyle(obj.first, (plot.show_ratio) ? 0.6666 : 1.);
      setAxisTitles(obj.first, plot);
    }

    gPad->Modified();
    gPad->Update();

    // We have the plot range. Compute the shaded area corresponding to the blinded area, if any
    if (!CommandLineCfg::get().unblind && h_data.get() && plot.blinded_range.valid()) {
        float x_start = plot.blinded_range.start;
        float x_end = plot.blinded_range.end;

        float y_start = gPad->GetUymin();
        float y_end = gPad->GetUymax();

        std::string options = "NB";

        if (plot.log_y) {
            options = options + " NDC";

            float lm = gPad->GetLeftMargin();
            float rm = 1. - gPad->GetRightMargin();
            float tm = 1. - gPad->GetTopMargin();
            float bm = gPad->GetBottomMargin();

            if (plot.log_x) {
                Range x_range = getXRange(toDraw[0].first);

                x_start = (rm - lm) * ((std::log(x_start) - std::log(x_range.start)) / (std::log(x_range.end) - std::log(x_range.start))) + lm;
                x_end = (rm - lm) * ((std::log(x_end) - std::log(x_range.start)) / (std::log(x_range.end) - std::log(x_range.start))) + lm;
            } else {
                x_start = (rm - lm) * ((x_start - gPad->GetUxmin()) / (gPad->GetUxmax() - gPad->GetUxmin())) + lm;
                x_end = (rm - lm) * ((x_end - gPad->GetUxmin()) / (gPad->GetUxmax() - gPad->GetUxmin())) + lm;
            }

            y_start = bm;
            y_end = tm;
        }

        std::shared_ptr<TPave> blinded_area(new TPave(x_start, y_start, x_end, y_end, 0, options.c_str()));
        blinded_area->SetFillStyle(m_plotIt.getConfiguration().blinded_range_fill_style);
        blinded_area->SetFillColor(m_plotIt.getConfiguration().blinded_range_fill_color);

        TemporaryPool::get().add(blinded_area);
        blinded_area->Draw("same");
    }

    // Draw all the requested lines for this plot
    auto resolveLine = [&](Line& line) {
        Range x_range = getXRange(toDraw[0].first);

        float y_range_start = gPad->GetUymin();
        float y_range_end = gPad->GetUymax();

        if (std::isnan(line.start.x))
            line.start.x = x_range.start;

        if (std::isnan(line.start.y))
            line.start.y = y_range_start;

        if (std::isnan(line.end.x))
            line.end.x = x_range.end;

        if (std::isnan(line.end.y))
            line.end.y = y_range_end;
    };

    for (Line& line: plot.lines) {
      resolveLine(line);

      std::shared_ptr<TLine> l(new TLine(line.start.x, line.start.y, line.end.x, line.end.y));
      TemporaryPool::get().add(l);

      l->SetLineColor(line.style->line_color);
      l->SetLineWidth(line.style->line_width);
      l->SetLineStyle(line.style->line_type);

      l->Draw("same");
    }

    // Redraw only axis
    toDraw[0].first->Draw("axis same");

    if (plot.show_ratio) {

      // Compute ratio and draw it
      low_pad->cd();
      low_pad->SetGridy();

      std::shared_ptr<TH1> h_low_pad_axis(static_cast<TH1*>(h_data->Clone()));
      h_low_pad_axis->SetDirectory(nullptr);
      h_low_pad_axis->Reset(); // Keep binning
      setRange(h_low_pad_axis.get(), plot.x_axis_range, plot.ratio_y_axis_range);

      setDefaultStyle(h_low_pad_axis.get(), 1. / 0.3333);
      h_low_pad_axis->GetYaxis()->SetTickLength(0.04);
      h_low_pad_axis->GetYaxis()->SetNdivisions(505, true);
      h_low_pad_axis->GetXaxis()->SetTickLength(0.07);

      h_low_pad_axis->Draw();

      std::shared_ptr<TGraphAsymmErrors> ratio = getRatio(h_data.get(), mc_histo_stat_only.get());
      ratio->Draw("P0 same");

      // Compute systematic errors in %
      std::shared_ptr<TH1> h_systematics(static_cast<TH1*>(h_low_pad_axis->Clone()));
      h_systematics->SetDirectory(nullptr);
      h_systematics->Reset(); // Keep binning
      h_systematics->SetMarkerSize(0);

      bool has_syst = false;
      for (uint32_t i = 1; i <= (uint32_t) h_systematics->GetNbinsX(); i++) {

        if (mc_histo_syst_only->GetBinContent(i) == 0 || mc_histo_syst_only->GetBinError(i) == 0)
          continue;

        // relative error, delta X / X
        float syst = mc_histo_syst_only->GetBinError(i) / mc_histo_syst_only->GetBinContent(i);

        h_systematics->SetBinContent(i, 1);
        h_systematics->SetBinError(i, syst);

        has_syst = true;
      }

      if (has_syst) {
        h_systematics->SetFillStyle(m_plotIt.getConfiguration().error_fill_style);
        h_systematics->SetFillColor(m_plotIt.getConfiguration().error_fill_color);
        setRange(h_systematics.get(), plot.x_axis_range, {});
        h_systematics->Draw("E2");
      }

      h_low_pad_axis->Draw("same");

      if (plot.fit_ratio) {
        float xMin, xMax;
        if (plot.ratio_fit_range.valid()) {
          xMin = plot.ratio_fit_range.start;
          xMax = plot.ratio_fit_range.end;
        } else {
          xMin = h_low_pad_axis->GetXaxis()->GetBinLowEdge(1);
          xMax = h_low_pad_axis->GetXaxis()->GetBinUpEdge(h_low_pad_axis->GetXaxis()->GetLast());
        }

        std::shared_ptr<TF1> fct = std::make_shared<TF1>("fit_function", plot.ratio_fit_function.c_str(), xMin, xMax);
        fct->SetNpx(m_plotIt.getConfiguration().ratio_fit_n_points);

        TFitResultPtr fit_result = ratio->Fit(fct.get(), "SMRNEQ");
        if (fit_result->IsValid()) {
          std::shared_ptr<TH1> errors = std::make_shared<TH1D>("errors", "errors", m_plotIt.getConfiguration().ratio_fit_n_points, xMin, xMax);
          errors->SetDirectory(nullptr);
          (TVirtualFitter::GetFitter())->GetConfidenceIntervals(errors.get(), 0.68);
          errors->SetStats(false);
          errors->SetMarkerSize(0);
          errors->SetFillColor(m_plotIt.getConfiguration().ratio_fit_error_fill_color);
          errors->SetFillStyle(m_plotIt.getConfiguration().ratio_fit_error_fill_style);
          errors->Draw("e3 same");

          fct->SetLineWidth(m_plotIt.getConfiguration().ratio_fit_line_width);
          fct->SetLineColor(m_plotIt.getConfiguration().ratio_fit_line_color);
          fct->SetLineStyle(m_plotIt.getConfiguration().ratio_fit_line_style);
          fct->Draw("same");

          if (plot.ratio_fit_legend.length() > 0) {
            uint32_t fit_parameters = fct->GetNpar();
            boost::format formatter = get_formatter(plot.ratio_fit_legend);

            for (uint32_t i = 0; i < fit_parameters; i++) {
              formatter % fct->GetParameter(i);
            }

            std::string legend = formatter.str();

            std::shared_ptr<TLatex> t(new TLatex(plot.ratio_fit_legend_position.x, plot.ratio_fit_legend_position.y, legend.c_str()));
            t->SetNDC(true);
            t->SetTextFont(43);
            t->SetTextSize(LABEL_FONTSIZE - 4);
            t->Draw();

            TemporaryPool::get().add(t);
          }

          TemporaryPool::get().add(errors);
          TemporaryPool::get().add(fct);
        }
      }

      h_low_pad_axis->Draw("same");
      ratio->Draw("P0 same");

      // Hide top pad label
      hideXTitle(toDraw[0].first);

      TemporaryPool::get().add(h_low_pad_axis);
      TemporaryPool::get().add(ratio);
      TemporaryPool::get().add(h_systematics);
      TemporaryPool::get().add(hi_pad);
      TemporaryPool::get().add(low_pad);
    }

    if (plot.fit) {
      float xMin, xMax;
      if (plot.fit_range.valid()) {
        xMin = plot.fit_range.start;
        xMax = plot.fit_range.end;
      } else {
        xMin = mc_stack->GetXaxis()->GetBinLowEdge(1);
        xMax = mc_stack->GetXaxis()->GetBinUpEdge(mc_stack->GetXaxis()->GetLast());
      }
      std::shared_ptr<TF1> fct = std::make_shared<TF1>("fit_function", plot.fit_function.c_str(), xMin, xMax);
      fct->SetNpx(m_plotIt.getConfiguration().fit_n_points);

      TH1* mc_hist = static_cast<TH1*>(mc_stack->GetStack()->At(mc_stack->GetNhists() - 1));
      TFitResultPtr fit_result = mc_hist->Fit(fct.get(), "SMRNEQ");
      if (fit_result->IsValid()) {
        std::shared_ptr<TH1> errors = std::make_shared<TH1D>("errors", "errors", m_plotIt.getConfiguration().fit_n_points, xMin, xMax);
        errors->SetDirectory(nullptr);
        (TVirtualFitter::GetFitter())->GetConfidenceIntervals(errors.get(), 0.68);
        errors->SetStats(false);
        errors->SetMarkerSize(0);
        errors->SetFillColor(m_plotIt.getConfiguration().fit_error_fill_color);
        errors->SetFillStyle(m_plotIt.getConfiguration().fit_error_fill_style);
        errors->Draw("e3 same");

        fct->SetLineWidth(m_plotIt.getConfiguration().fit_line_width);
        fct->SetLineColor(m_plotIt.getConfiguration().fit_line_color);
        fct->SetLineStyle(m_plotIt.getConfiguration().fit_line_style);
        fct->Draw("same");

        if (plot.fit_legend.length() > 0) {
          uint32_t fit_parameters = fct->GetNpar();
          boost::format formatter = get_formatter(plot.fit_legend);

          for (uint32_t i = 0; i < fit_parameters; i++) {
            formatter % fct->GetParameter(i);
          }

          std::string legend = formatter.str();

          std::shared_ptr<TLatex> t(new TLatex(plot.fit_legend_position.x, plot.fit_legend_position.y, legend.c_str()));
          t->SetNDC(true);
          t->SetTextFont(43);
          t->SetTextSize(LABEL_FONTSIZE - 4);
          t->Draw();

          TemporaryPool::get().add(t);
        }

        TemporaryPool::get().add(errors);
        TemporaryPool::get().add(fct);
      }
    }

    gPad->Modified();
    gPad->Update();
    gPad->RedrawAxis();

    if (hi_pad.get())
      hi_pad->cd();

    return global_summary;
  }

  void TH1Plotter::setHistogramStyle(const File& file) {
    TH1* h = dynamic_cast<TH1*>(file.object);

    std::shared_ptr<PlotStyle> style = m_plotIt.getPlotStyle(file);

    if (style->fill_color != -1)
      h->SetFillColor(style->fill_color);

    if (style->fill_type != -1)
      h->SetFillStyle(style->fill_type);

    if (style->line_color != -1)
      h->SetLineColor(style->line_color);

    if (style->line_width != -1)
      h->SetLineWidth(style->line_width);

    if (style->line_type != -1)
      h->SetLineStyle(style->line_type);

    if (style->marker_size != -1)
      h->SetMarkerSize(style->marker_size);

    if (style->marker_color != -1)
      h->SetMarkerColor(style->marker_color);

    if (style->marker_type != -1)
      h->SetMarkerStyle(style->marker_type);

    if (file.type == MC && style->line_color == -1 && style->fill_color != -1)
      h->SetLineColor(style->fill_color);
  }

  void TH1Plotter::addOverflow(TH1* h, Type type, const Plot& plot) {

    if (!h || !h->GetEntries())
        return;

    size_t first_bin = 1;
    size_t last_bin = h->GetNbinsX();

    if (plot.x_axis_range.valid()) {
      std::shared_ptr<TH1> copy(dynamic_cast<TH1*>(h->Clone()));
      copy->SetDirectory(nullptr);
      copy->GetXaxis()->SetRangeUser(plot.x_axis_range.start, plot.x_axis_range.end);

      // Find first and last bin corresponding to the given range
      first_bin = copy->GetXaxis()->GetFirst();
      last_bin = copy->GetXaxis()->GetLast();
    }

    // GetBinError returns sqrt(SumW2) for a given bin
    // SetBinError updates SumW2 for a given bin with error*error

    float underflow = 0;
    float underflow_sumw2 = 0;
    for (size_t i = 0; i < first_bin; i++) {
      underflow += h->GetBinContent(i);
      underflow_sumw2 += (h->GetBinError(i) * h->GetBinError(i));
    }
    float overflow = 0;
    float overflow_sumw2 = 0;
    for (size_t i = last_bin + 1; i <= (size_t) h->GetNbinsX() + 1; i++) {
      overflow += h->GetBinContent(i);
      overflow_sumw2 += (h->GetBinError(i) * h->GetBinError(i));
    }
    // Clear out-of-range bin content so that Integral() still returns the right value
    for (size_t i = 1; i < first_bin; i++) {
      h->SetBinContent(i, 0);
    }
    for (size_t i = last_bin + 1; i < (size_t) h->GetNbinsX() + 1; i++) {
      h->SetBinContent(i, 0);
    }
    // Clear also underflow and overflow bins (SetBinContent on these may try to extend the axes)
    h->ClearUnderflowAndOverflow();

    float first_bin_content = h->GetBinContent(first_bin);
    float first_bin_sumw2 = h->GetBinError(first_bin) * h->GetBinError(first_bin);

    float last_bin_content = h->GetBinContent(last_bin);
    float last_bin_sumw2 = h->GetBinError(last_bin) * h->GetBinError(last_bin);

    h->SetBinContent(first_bin, first_bin_content + underflow);
    if (type != DATA)
        h->SetBinError(first_bin, sqrt(underflow_sumw2 + first_bin_sumw2));

    h->SetBinContent(last_bin, last_bin_content + overflow);
    if (type != DATA)
        h->SetBinError(last_bin, sqrt(overflow_sumw2 + last_bin_sumw2));
  }
}
