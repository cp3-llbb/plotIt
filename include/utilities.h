#pragma once

#include <defines.h>

#include <plotIt.h>

#include <boost/format.hpp>

namespace plotIt {
  TStyle* createStyle();

  boost::format get_formatter(const std::string format_string);

  #define CAST_AND_CALL(OBJECT, FUNCTION, ...) \
      if (dynamic_cast<TH1*>(OBJECT)) \
        FUNCTION(dynamic_cast<TH1*>(OBJECT), ##__VA_ARGS__); \
      else if (dynamic_cast<THStack*>(OBJECT)) \
        FUNCTION(dynamic_cast<THStack*>(OBJECT), ##__VA_ARGS__);

  #define CAST_AND_RETURN(OBJECT, FUNCTION, ...) \
      if (dynamic_cast<TH1*>(OBJECT)) \
        return FUNCTION(dynamic_cast<TH1*>(OBJECT), ##__VA_ARGS__); \
      else if (dynamic_cast<THStack*>(OBJECT)) \
        return FUNCTION(dynamic_cast<THStack*>(OBJECT), ##__VA_ARGS__);

  #define CAST_TO_HIST_AND_CALL(OBJECT, FUNCTION, ...) \
      if (dynamic_cast<TH1*>(OBJECT)) \
        FUNCTION(dynamic_cast<TH1*>(OBJECT), ##__VA_ARGS__); \
      else if (dynamic_cast<THStack*>(OBJECT)) \
        FUNCTION(dynamic_cast<THStack*>(OBJECT)->GetHistogram(), ##__VA_ARGS__);

  #define CAST_TO_HIST_AND_RETURN(OBJECT, FUNCTION, ...) \
      if (dynamic_cast<TH1*>(OBJECT)) \
        return FUNCTION(dynamic_cast<TH1*>(OBJECT), ##__VA_ARGS__); \
      else if (dynamic_cast<THStack*>(OBJECT)) \
        return FUNCTION(dynamic_cast<THStack*>(OBJECT)->GetHistogram(), ##__VA_ARGS__);

  #define ADD_PAIRS(PAIR1, PAIR2) \
      PAIR1.first += PAIR2.first; PAIR1.second += PAIR2.second;

  template<class T>
    void setAxisTitles(T* object, Plot& plot) {
      if (plot.x_axis.length() > 0 && object->GetXaxis()) {
        object->GetXaxis()->SetTitle(plot.x_axis.c_str());
      }

      if (plot.y_axis.length() > 0 && object->GetYaxis()) {
        float binSize = object->GetXaxis()->GetBinWidth(1);
        std::string title = plot.y_axis;

        boost::format formatter = get_formatter(plot.y_axis_format);
        object->GetYaxis()->SetTitle((formatter % title % binSize).str().c_str());
      }

      if (plot.show_ratio && object->GetXaxis())
        object->GetXaxis()->SetLabelSize(0);
    }

  void setAxisTitles(TObject* object, Plot& plot);

  template<class T>
    void setDefaultStyle(T* object, float topBottomScaleFactor) {

      // Remove title
      object->SetBit(TH1::kNoTitle);
      object->SetLabelFont(43, "XYZ");
      object->SetTitleFont(43, "XYZ");
      object->SetLabelSize(LABEL_FONTSIZE, "XYZ");
      object->SetTitleSize(TITLE_FONTSIZE, "XYZ");

      object->GetYaxis()->SetTitle("Data / MC");
      object->GetYaxis()->SetNdivisions(510);
      object->GetYaxis()->SetTitleOffset(2.5);
      object->GetYaxis()->SetLabelOffset(0.01);
      object->GetYaxis()->SetTickLength(0.03);

      object->GetXaxis()->SetTitleOffset(1.5 * topBottomScaleFactor);
      object->GetXaxis()->SetLabelOffset(0.012 * topBottomScaleFactor);
      object->GetXaxis()->SetTickLength(0.03);

      // No stats box
      object->SetStats(false);
      
    }

  void setDefaultStyle(TObject* object, float topBottomScaleFactor);


  template<class T>
    void hideXTitle(T* object) {
      object->GetXaxis()->SetTitle("");
      object->GetXaxis()->SetTitleSize();
    }

  void hideXTitle(TObject* object);

  template<class T>
    float getMaximum(T* object) {
      return object->GetMaximum();
    }

  float getMaximum(TObject* object);

  template<class T>
    float getMinimum(T* object) {
      return object->GetMinimum();
    }

  float getMinimum(TObject* object);

  template<class T>
    void setMaximum(T* object, float minimum) {
      object->SetMaximum(minimum);
    }

  void setMaximum(TObject* object, float minimum);

  template<class T>
    void setMinimum(T* object, float minimum) {
      object->SetMinimum(minimum);
    }

  void setMinimum(TObject* object, float minimum);

  template<class T>
    void setRange(T* object, const Range& x_range, const Range& y_range) {
      if (x_range.valid())
        object->GetXaxis()->SetRangeUser(x_range.start, x_range.end);
      if (y_range.valid()) {
        object->SetMinimum(y_range.start);
        object->SetMaximum(y_range.end);
      }
    }

  void setRange(TObject* object, const Range& x_range, const Range& y_range);

  template<class T>
    Range getXRange(T* object) {
      Range range;
      range.start = object->GetXaxis()->GetBinLowEdge(object->GetXaxis()->GetFirst());
      range.end = object->GetXaxis()->GetBinUpEdge(object->GetXaxis()->GetLast());

      return range;
    }

  Range getXRange(TObject* object);

  template<class T>
  float getPositiveMinimum(T* object) {
    return object->GetMinimum(0);
  }

  float getPositiveMinimum(TObject* object);

  // replace all occurences of "old" in "s" by "rep"
  void replace_substr(std::string &s, const std::string &old, const std::string &rep);

  int16_t loadColor(const YAML::Node& node);

  inline std::vector<std::string> glob(const std::string& pat) {
      glob_t glob_result;
      glob(pat.c_str(), GLOB_TILDE, NULL, &glob_result);

      std::vector<std::string> ret;
      for(unsigned int i = 0;i < glob_result.gl_pathc; ++i){
          ret.push_back(std::string(glob_result.gl_pathv[i]));
      }

      globfree(&glob_result);
      return ret;
  }

  inline std::string truncate(const std::string& str, size_t max_len) {
      if (str.length() > max_len - 1) {
          std::string ret = str;
          ret.resize(max_len - 1);
          return ret + u8"…";
      } else {
          return str;
      }
  }

  TDirectory* getDirectory(TDirectoryFile* root, const boost::filesystem::path& directory, bool create = true);

    std::string applyRenaming(const std::vector<RenameOp>& ops, const std::string input);
}
