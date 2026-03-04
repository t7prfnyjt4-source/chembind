import React, { useMemo } from "react";
import type { JobDoc, JobStatus } from "../hooks/useBatchJob";

function statusStyle(status: JobStatus) {
  switch (status) {
    case "queued":
      return { bg: "#eef2ff", fg: "#3730a3", label: "Queued" };
    case "running":
      return { bg: "#ecfeff", fg: "#155e75", label: "Running" };
    case "completed":
      return { bg: "#ecfdf5", fg: "#065f46", label: "Completed" };
    case "failed":
      return { bg: "#fef2f2", fg: "#991b1b", label: "Failed" };
  }
}

export function BatchJobProgress({ job }: { job: JobDoc & { id: string } }) {
  const s = statusStyle(job.status);

  const total = job.total ?? 0;
  const processed = job.processed ?? 0;

  const pct = useMemo(() => {
    if (!total || total <= 0) return job.status === "completed" ? 100 : 0;
    const raw = Math.round((processed / total) * 100);
    return Math.max(0, Math.min(100, raw));
  }, [total, processed, job.status]);

  const success = job.successCount ?? 0;
  const failure = job.failureCount ?? 0;

  return (
    <div style={{ border: "1px solid #e5e7eb", borderRadius: 12, padding: 12, marginTop: 12 }}>
      <div style={{ display: "flex", justifyContent: "space-between", gap: 12, alignItems: "center" }}>
        <div style={{ display: "flex", alignItems: "center", gap: 10 }}>
          <span style={{ padding: "4px 10px", borderRadius: 999, background: s.bg, color: s.fg, fontWeight: 600 }}>
            {s.label}
          </span>
          <span style={{ color: "#6b7280" }}>
            Job: <code>{job.id}</code>
          </span>
        </div>

        <div style={{ color: "#374151", fontVariantNumeric: "tabular-nums" }}>
          {processed}/{total || "?"} processed • ✅ {success} • ❌ {failure}
        </div>
      </div>

      <div style={{ marginTop: 10 }}>
        <div style={{ height: 10, background: "#f3f4f6", borderRadius: 999, overflow: "hidden" }}>
          <div style={{ width: `${pct}%`, height: "100%", background: "#111827" }} />
        </div>
        <div style={{ marginTop: 6, color: "#6b7280" }}>{pct}%</div>
      </div>

      {job.status === "failed" && job.errorMessage && (
        <div style={{ marginTop: 10, padding: 10, borderRadius: 10, background: "#fef2f2", color: "#991b1b" }}>
          <strong>Error:</strong> {job.errorMessage}
        </div>
      )}
    </div>
  );
}
