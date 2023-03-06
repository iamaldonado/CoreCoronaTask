# How to run CoreCorona task

This working example gets a toy hyperon Global Polarization. It is based on the lego trains introduced by V.Riabov [presentation](https://indico.jinr.ru/event/3391/contributions/18477/attachments/13910/23277/AnalysisFramework_RiabovV.pdf) at MPD Cross-PWG Meeting.

It employs the classes MpdAnalysisManager, MpdAnalysisEvent, and MpdAnalysisTask commited by D. Peresunko originally.

## Task Basic Structure

Tasks called by the MpdAnalysis Manager should be derived from MpdAnalysisTask and implements the following methods 

- ** Output
Users should perapre objects to fill (histograms, trees, etc.)
```ruby
	void UserInit();
```
- ** Main
Method is called for each event, data are provided by container MpdAnalysisEvent
```ruby
	void ProcessEvent(MpdAnalysisEvent &event);
```
- ** End
Method is called when scan is finished but class data are not written yet
```ruby
	void Finish();
```
## MpdAnalysisEvent Class

Contains references to all branched in the DST file for this event and may contain extra global variables of interest (centrality, event plane, T_0, n-sigma matching ... etc.)





 

