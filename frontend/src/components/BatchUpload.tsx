import React, { useMemo, useState } from "react";
import { getAuth } from "firebase/auth";

import { parseSmilesCsv } from "../utils/csvParser";
import { sha256Hex } from "../utils/sha256";

import { useBatchJob } from "../hooks/useBatchJob";
import { useJobList } from "../hooks/useJobList";
import { useJobItems } from "../hooks/useJobItems";
import { BatchJobProgress } from "./BatchJobProgress";
import { JobResultsTable } from "./JobResultsTable";

type Props = {
  token: string;
};

export default function BatchUpload({ token }: Props) {
  const auth = getAuth();
  const uid = auth.currentUser?.uid ?? null;

  const [csvText, setCsvText] = useState("");
  const [rows, setRows] = useState<any[]>([]);
  const [jobId, setJobId] = useState<string | null>(null);
  const [viewJobId, setViewJobId] = useState<string | null>(null);
  const [error, setError] = useState<string | null>(null);

  const preview = useMemo(() => rows.slice(0, 10), [rows]);

  const { job: activeJob, loading: activeJobLoading, error: activeJobError } = useBatchJob(uid, jobId);
  const { jobs, loading: jobsLoading, error: jobsError } = useJobList(uid);
  const { items, loading: itemsLoading, error: itemsError } = useJobItems(uid, viewJobId, Boolean(viewJobId));

  async function onFile(file: File) {
    setError(null);
    setJobId(null);
    setViewJobId(null);

    const text = await file.text();
    setCsvText(text);

    const parsed = parseSmilesCsv(text);
    setRows(parsed);

    if (parsed.length === 0) {
      setError("No SMILES rows found.");
    }
  }

  async function submit() {
    if (!token || rows.length === 0) return;

    try {
      setError(null);

      const key = await sha256Hex(csvText);

      const res = await fetch("/api/batch", {
        method: "POST",
        headers: {
          "Content-Type": "application/json",
          Authorization: `Bearer ${token}`,
          "Idempotency-Key": key,
        },
        body: JSON.stringify({
          rows,
          source: "csv",
        }),
      });

      if (!res.ok) {
        const msg = await res.text().catch(() => "");
        throw new Error(msg || `HTTP ${res.status}`);
      }

      const data = await res.json();
      const newJobId = data.job_id as string | undefined;

      if (!newJobId) {
        throw new Error("Missing job_id in response.");
      }

      setJobId(newJobId);     // track progress
      setViewJobId(null);     // close old results
    } catch (e) {
      setError("Upload failed.");
    }
  }

  return (
    <div>
      <h2>Batch Upload</h2>

      {!uid && (
        <p style={{ color: "#6b7280" }}>
          Sign in to see job progress/results.
        </p>
      )}

      <input
        aria-label="csv-file"
        type="file"
        accept=".csv"
        onChange={(e) => {
          const f = e.target.files?.[0];
          if (f) onFile(f);
        }}
      />

      {rows.length > 0 && (
        <div aria-label="preview">
          <p>
            Preview (first {Math.min(10, rows.length)} of {rows.length})
          </p>
          <table>
            <thead>
              <tr>
                <th>rowId</th>
                <th>smiles</th>
              </tr>
            </thead>
            <tbody>
              {preview.map((r) => (
                <tr key={r.rowId}>
                  <td>{r.rowId}</td>
                  <td>{r.smiles}</td>
                </tr>
              ))}
            </tbody>
          </table>
        </div>
      )}

      <button aria-label="submit" onClick={submit} disabled={rows.length === 0}>
        Submit
      </button>

      {jobId && (
        <p aria-label="job-id">
          job_id: {jobId}
        </p>
      )}

      {activeJobLoading && <div style={{ color: "#6b7280", marginTop: 8 }}>Loading job…</div>}
      {activeJobError && <div style={{ color: "#991b1b", marginTop: 8 }}>{activeJobError}</div>}

      {activeJob && (
        <>
          <BatchJobProgress job={activeJob} />

          {activeJob.status === "completed" && (
            <div style={{ marginTop: 10 }}>
              <button
                onClick={() => setViewJobId(activeJob.id)}
                style={{
                  padding: "8px 12px",
                  borderRadius: 10,
                  border: "1px solid #e5e7eb",
                  background: "#fff",
                }}
              >
                View Results 📄
              </button>
            </div>
          )}
        </>
      )}

      <div style={{ marginTop: 18 }}>
        <h3 style={{ margin: 0 }}>Past jobs</h3>

        {jobsLoading && <div style={{ color: "#6b7280", marginTop: 6 }}>Loading jobs…</div>}
        {jobsError && <div style={{ color: "#991b1b", marginTop: 6 }}>{jobsError}</div>}

        {!jobsLoading && jobs.length === 0 && (
          <div style={{ color: "#6b7280", marginTop: 6 }}>
            No jobs yet. Upload a CSV to start 🚀
          </div>
        )}

        <div style={{ display: "grid", gap: 8, marginTop: 10 }}>
          {jobs.map((j) => (
            <div
              key={j.id}
              style={{
                border: "1px solid #e5e7eb",
                borderRadius: 12,
                padding: 10,
                display: "flex",
                justifyContent: "space-between",
                gap: 10,
              }}
            >
              <div>
                <div style={{ fontWeight: 600 }}>
                  <code>{j.id}</code>
                </div>
                <div style={{ color: "#6b7280", fontSize: 13 }}>
                  {j.status} • {j.processed ?? 0}/{j.total ?? "?"}
                </div>
              </div>

              <div style={{ display: "flex", gap: 8, alignItems: "center" }}>
                <button
                  onClick={() => {
                    setJobId(j.id);
                    setViewJobId(null);
                  }}
                  style={{
                    padding: "6px 10px",
                    borderRadius: 10,
                    border: "1px solid #e5e7eb",
                    background: "#fff",
                  }}
                >
                  Track
                </button>

                <button
                  disabled={j.status !== "completed"}
                  onClick={() => setViewJobId(j.id)}
                  style={{
                    padding: "6px 10px",
                    borderRadius: 10,
                    border: "1px solid #e5e7eb",
                    background: j.status === "completed" ? "#fff" : "#f9fafb",
                    color: j.status === "completed" ? "#111827" : "#9ca3af",
                    cursor: j.status === "completed" ? "pointer" : "not-allowed",
                  }}
                >
                  View Results
                </button>
              </div>
            </div>
          ))}
        </div>
      </div>

      {viewJobId && (
        <div style={{ marginTop: 18 }}>
          <h3 style={{ margin: 0 }}>
            Results for <code>{viewJobId}</code>
          </h3>

          {itemsLoading && <div style={{ color: "#6b7280", marginTop: 6 }}>Loading results…</div>}
          {itemsError && <div style={{ color: "#991b1b", marginTop: 6 }}>{itemsError}</div>}
          {!itemsLoading && !itemsError && <JobResultsTable items={items} />}
        </div>
      )}

      {error && (
        <p aria-label="error" style={{ color: "red" }}>
          {error}
        </p>
      )}
    </div>
  );
}
