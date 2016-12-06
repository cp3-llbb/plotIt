#pragma once

#include <boost/algorithm/string/join.hpp>
#include <boost/filesystem.hpp>
#include <boost/optional.hpp>

#include <iostream>

#include <defines.h>
#include <uuid.h>
#include <systematics.h>

#include <yaml-cpp/yaml.h>

#include <TObject.h>
#include <TFile.h>
#include <TChain.h>

class TLegendEntry;

namespace plotIt {

  enum Type {
    MC,
    SIGNAL,
    DATA
  };

  enum Location {
    TOP,
    BOTTOM
  };

  inline Type string_to_type(const std::string& type) {
      if (type == "signal")
          return SIGNAL;
      else if (type == "data")
          return DATA;
      else
          return MC;
  }

  inline std::string type_to_string(const Type& type) {
      switch (type) {
        case MC:
            return "Simulation";

        case SIGNAL:
            return "Signal";

        case DATA:
            return "Data";
      }

      return "Unknown";
  }

  enum ErrorsType {
      Normal = 0,
      Poisson = 1,
      Poisson2 = 2
  };

  inline ErrorsType string_to_errors_type(const std::string& s) {
      ErrorsType errors_type;

      if (s == "normal")
          errors_type = Normal;
      else if (s == "poisson2")
          errors_type = Poisson2;
      else 
          errors_type = Poisson;

      return errors_type;
  }

  enum Log {
    False = 0,
    True,
    Both
  };

  inline Log parse_log(const YAML::Node& node) {
    if (node.as<std::string>() == "both")
      return Both;
    else
      return node.as<bool>() ? True : False;
  }

  enum Orientation {
      UNSPECIFIED,
      HORIZONTAL,
      VERTICAL
  };

  struct LineStyle {
    float line_width = 1;
    int16_t line_color = 1;
    int16_t line_type = 1;

    LineStyle() = default;
    LineStyle(const YAML::Node& node);

    void parse(const YAML::Node& node);
  };

  struct PlotStyle: public LineStyle {

    // Style
    float marker_size;
    int16_t marker_color;
    int16_t marker_type;
    int16_t fill_color;
    int16_t fill_type;
    std::string drawing_options;

    // Legend
    std::string legend;
    std::string legend_style;
    int16_t legend_order = 0;

    void loadFromYAML(const YAML::Node& node, Type type);
  };

  struct RenameOp {
      std::regex from;
      std::string to;
  };

  struct File {
    std::string path;
    std::string pretty_name;
    size_t id = 0;

    // For MC and Signal
    float cross_section = 1.;
    float branching_ratio = 1.;
    float generated_events = 1.;
    float scale = 1.;

    // Only MC files with the same stack index will be
    // merged together
    int64_t stack_index = 0;

    std::shared_ptr<PlotStyle> plot_style;
    std::string legend_group;
    std::string yields_group;

    Type type = MC;

    TObject* object = nullptr;
    std::map<std::string, TObject*> objects;

    std::vector<SystematicSet>* systematics;
    std::map<std::string, std::vector<SystematicSet>> systematics_cache;

    int16_t order = std::numeric_limits<int16_t>::min();

    std::shared_ptr<TChain> chain;

    std::shared_ptr<TFile> handle;
    std::map<std::string, std::shared_ptr<TFile>> friend_handles;

    // Renaming
    std::vector<RenameOp> renaming_ops;
  };

  struct Group {
    std::string name;
    std::shared_ptr<PlotStyle> plot_style;

    bool added = false;
  };

  struct Point {
    float x = std::numeric_limits<float>::quiet_NaN();
    float y = std::numeric_limits<float>::quiet_NaN();

    bool operator==(const Point& other) {
      return
        (fabs(x - other.x) < 1e-6) &&
        (fabs(y - other.y) < 1e-6);
    }

    bool valid() const {
      return !std::isnan(x) && !std::isnan(y);
    }

    Point() = default;
    Point(std::initializer_list<float> c) {
      assert(c.size() == 2);
      x = *c.begin();
      y = *(c.begin() + 1);
    }
  };

  struct Range {
    float start = std::numeric_limits<float>::quiet_NaN();
    float end = std::numeric_limits<float>::quiet_NaN();

    bool operator==(const Range& other) {
      return
        (std::abs(start - other.start) < 1e-6) &&
        (std::abs(end - other.end) < 1e-6);
    }

    bool valid() const {
      return !std::isnan(start) && !std::isnan(end);
    }

    Range() = default;
    Range(std::initializer_list<float> c) {
      assert(c.size() == 2);
      start = *c.begin();
      end = *(c.begin() + 1);
    }
  };

  struct Position {
    float x1 = 0;
    float y1 = 0;

    float x2 = 0;
    float y2 = 0;

    Position(float x1, float y1, float x2, float y2):
        x1(x1), y1(y1), x2(x2), y2(y2) {
        // Empty
    }

    Position() = default;

    bool operator==(const Position& other) {
      return
        (fabs(x1 - other.x1) < 1e-6) &&
        (fabs(y1 - other.y1) < 1e-6) &&
        (fabs(x2 - other.x2) < 1e-6) &&
        (fabs(y2 - other.y2) < 1e-6);
    }
  };

  struct Label {
    std::string text;
    uint32_t size = LABEL_FONTSIZE;
    Point position;
  };

  struct Line {
    Point start;
    Point end;

    boost::optional<LineStyle> style;
    Location pad = TOP;

    bool operator==(const Line& other) {
      return ((start == other.start) && (end == other.end));
    }

    bool valid() const {
      return start.valid() && end.valid();
    }

    Line() = default;
    Line(std::initializer_list<Point> c) {
      assert(c.size() == 2);
      start = *c.begin();
      end = *(c.begin() + 1);
    }

    Line(const YAML::Node& node, Orientation);
  };

  struct Plot {
    std::string name;
    std::string output_suffix;
    std::string uid = get_uuid();
    std::string exclude;
    std::string book_keeping_folder;
    std::vector<RenameOp> renaming_ops;

    bool no_data = false;
    bool override = false; // flag to plot only those which have it true (if at least one plot has it true)
    bool normalized = false;
    bool log_y = false;
    bool log_x = false;

    std::string x_axis;
    std::string y_axis = "Events";
    std::string y_axis_format;
    bool y_axis_show_zero = false;
    std::string ratio_y_axis_title = "Data / MC";

    // Axis range
    Range x_axis_range;
    Range log_x_axis_range;
    Range y_axis_range;
    Range log_y_axis_range;
    Range ratio_y_axis_range = {0.5, 1.5};

    // Blind range
    Range blinded_range;

    uint16_t binning_x;  // Only used in tree mode
    uint16_t binning_y;  // Only used in tree mode

    std::string draw_string;  // Only used in tree mode
    std::string selection_string;  // Only used in tree mode

    std::vector<std::string> save_extensions = {"pdf"};

    bool show_ratio = false;

    bool fit = false;
    std::string fit_function = "gaus";
    std::string fit_legend = "#scale[1.6]{#splitline{#mu = %2$.3f}{#sigma = %3$.3f}}";
    Point fit_legend_position = {0.22, 0.87};
    Range fit_range;

    bool fit_ratio = false;
    std::string ratio_fit_function = "pol1";
    std::string ratio_fit_legend;
    Point ratio_fit_legend_position = {0.20, 0.38};
    Range ratio_fit_range;

    bool show_errors = true;
    bool show_overflow = false;

    std::string inherits_from = "TH1";

    uint16_t rebin = 1;

    std::vector<Label> labels;

    std::string extra_label;

    Position legend_position;
    size_t legend_columns;

    ErrorsType errors_type = Poisson;

    bool use_for_yields = false;
    std::string yields_title;
    int yields_table_order = 0;

    bool is_rescaled = false;

    bool sort_by_yields = true;

    std::vector<Line> lines;

    // Axis label size
    float x_axis_label_size = LABEL_FONTSIZE;
    float y_axis_label_size = LABEL_FONTSIZE;

    // Show or hide ticks for each axis
    bool x_axis_hide_ticks = false;
    bool y_axis_hide_ticks = false;
    
    void print() {
      std::cout << "Plot '" << name << "'" << std::endl;
      std::cout << "\tx_axis: " << x_axis << std::endl;
      std::cout << "\ty_axis: " << y_axis << std::endl;
      std::cout << "\tshow_ratio: " << show_ratio << std::endl;
      std::cout << "\tinherits_from: " << inherits_from << std::endl;
      std::cout << "\tsave_extensions: " << boost::algorithm::join(save_extensions, ", ") << std::endl;
    }

    Plot Clone(const std::string& new_name) {
      Plot clone = *this;
      clone.name = new_name;
      clone.uid = get_uuid();

      return clone;
    }
  };

  struct Legend {
    Position position = {0.6, 0.6, 0.9, 0.9};
    size_t columns = 1;
  };

  struct LegendEntry {
      TObject* object = nullptr;
      std::string legend;
      std::string style;
      int16_t order = 0;

      int16_t fill_style = 0;
      int16_t fill_color = 0;
      uint16_t line_width = 0;

      LegendEntry() = default;
      LegendEntry(TObject* object, const std::string& legend, const std::string& style, int16_t order);
      LegendEntry(const std::string& legend, const std::string& style, int16_t fill_style, int16_t fill_color, uint16_t line_width);
      void stylize(TLegendEntry* entry);
  };

  struct Configuration {
    float width = 800;
    float height = 800;
    float luminosity = -1;
    float scale = 1;
    bool no_lumi_rescaling = false;

    // Systematics
    float luminosity_error_percent = 0;

    std::string y_axis_format = "%1% / %2$.2f";
    std::string ratio_y_axis_title = "Data / MC";

    int16_t error_fill_color = 42;
    int16_t error_fill_style = 3154;

    uint16_t fit_n_points = 1000;
    int16_t fit_line_color = 46;
    int16_t fit_line_width = 1;
    int16_t fit_line_style = 1;
    int16_t fit_error_fill_color = 42;
    int16_t fit_error_fill_style = 1001;

    uint16_t ratio_fit_n_points = 1000;
    int16_t ratio_fit_line_color = 46;
    int16_t ratio_fit_line_width = 1;
    int16_t ratio_fit_line_style = 1;
    int16_t ratio_fit_error_fill_color = 42;
    int16_t ratio_fit_error_fill_style = 1001;

    LineStyle line_style;

    std::vector<Label> labels;

    std::string experiment = "CMS";
    std::string extra_label;

    std::string lumi_label;

    std::string root = "./";

    bool show_overflow = false;
    bool transparent_background = false;

    std::string mode = "hist"; // "tree" or "hist"
    std::string tree_name;

    ErrorsType errors_type = Poisson;

    float yields_table_stretch = 1.15;
    std::string yields_table_align = "h";
    std::string yields_table_text_align = "c";
    int yields_table_num_prec_yields = 1;
    int yields_table_num_prec_ratio = 2;

    int16_t blinded_range_fill_color = 42;
    int16_t blinded_range_fill_style = 1001;

    std::string uncertainty_label = "Uncertainty";
    std::map<Type, std::vector<LegendEntry>> static_legend_entries;

    std::string book_keeping_file_name;
    std::shared_ptr<TFile> book_keeping_file;

    // Axis label size
    float x_axis_label_size = LABEL_FONTSIZE;
    float y_axis_label_size = LABEL_FONTSIZE;

    // Show or not opposite axis ticks
    bool x_axis_top_ticks = true;
    bool y_axis_right_ticks = true;
  };
}

namespace YAML {
  template<>
    struct convert<plotIt::Position> {
      static Node encode(const plotIt::Position& rhs) {
        Node node;
        node.push_back(rhs.x1);
        node.push_back(rhs.y1);
        node.push_back(rhs.x1);
        node.push_back(rhs.y1);

        return node;
      }

      static bool decode(const Node& node, plotIt::Position& rhs) {
        if(!node.IsSequence() || node.size() != 4)
          return false;

        rhs.x1 = node[0].as<float>();
        rhs.y1 = node[1].as<float>();
        rhs.x2 = node[2].as<float>();
        rhs.y2 = node[3].as<float>();

        return true;
      }
    };

  template<>
    struct convert<plotIt::Point> {
      static Node encode(const plotIt::Point& rhs) {
        Node node;
        node.push_back(rhs.x);
        node.push_back(rhs.y);

        return node;
      }

      static bool decode(const Node& node, plotIt::Point& rhs) {
        if(!node.IsSequence() || node.size() != 2)
          return false;

        rhs.x = node[0].as<float>();
        rhs.y = node[1].as<float>();

        return true;
      }
    };

  template<>
    struct convert<plotIt::Range> {
      static Node encode(const plotIt::Range& rhs) {
        Node node;
        node.push_back(rhs.start);
        node.push_back(rhs.end);

        return node;
      }

      static bool decode(const Node& node, plotIt::Range& rhs) {
        if(!node.IsSequence() || node.size() != 2)
          return false;

        rhs.start = node[0].as<float>();
        rhs.end = node[1].as<float>();

        return true;
      }
    };
}

