import { useState, useEffect } from "react";
import { doc, onSnapshot } from "firebase/firestore";
import { db } from "../firebase/config";
import { auth } from "../firebase/config";

export interface DockingJob {
  jobId: string;
  status: string;
  smiles: string;
  center: number[];
  size: number[];
  exhaustiveness: number;
  numModes: number;
  poseCount?: number;
  bestScore?: number;
  error?: string;
  created_at: string;
  updated_at: string;
}

export function useDockingJob(jobId: string | null) {
  const [job, setJob] = useState<DockingJob | null>(null);

  useEffect(() => {
    const uid = auth.currentUser?.uid;
    if (!jobId || !uid) return;

    const docRef = doc(db, "users", uid, "docking_jobs", jobId);
    const unsub = onSnapshot(docRef, (snap) => {
      if (snap.exists()) {
        setJob(snap.data() as DockingJob);
      }
    });

    return () => unsub();
  }, [jobId]);

  return job;
}
