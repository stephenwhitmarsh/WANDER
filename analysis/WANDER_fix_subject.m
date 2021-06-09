% fix subject 17

isubject = 17;
force = 1;
rootpath = 1;
restingstate = 0;
timing = 'cue';
latency = [0 30];

WANDER_epoch_MEG(isubject,force,rootpath,restingstate)
WANDER_ICA(isubject,force,rootpath,restingstate)
WANDER_filter_EGG(isubject,force,rootpath,restingstate)
WANDER_behaviour(isubject,force,rootpath)
WANDER_TFR(isubject,force,'cue',rootpath,restingstate)
WANDER_artefact_detection_EGG(isubject,force,rootpath,restingstate)

WANDER_artefact_detection_MEG(isubject,force,rootpath,restingstate)
WANDER_blink_detection(isubject,force,rootpath,restingstate)
WANDER_ICA(isubject,force,rootpath,restingstate)
WANDER_artefact_detection_MEG_after_ICA(isubject,force,rootpath,restingstate)

WANDER_ERF(isubject,force,'cue',rootpath)
WANDER_TFR_phaselocking_balanced(isubject,force,timing,latency,rootpath)

WANDER_Heart(isubject,force,timing,rootpath)
