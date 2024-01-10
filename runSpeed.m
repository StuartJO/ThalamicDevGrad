function R = runSpeed(TotalDist,TotalTime,WalkDist,WalkSpeed)

R = (TotalDist-WalkDist) / ( (TotalTime/60) - (WalkDist/WalkSpeed) );