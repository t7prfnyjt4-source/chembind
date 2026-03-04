import React from "react";
import { describe, it, expect, vi, beforeEach } from "vitest";
import { render, screen } from "@testing-library/react";
import userEvent from "@testing-library/user-event";

import BatchUpload from "./BatchUpload";

// ---- Mocks ----

// Firebase auth mock (mutable so we can flip to "no user" in a test)
let mockAuthUser: { uid: string } | null = { uid: "testUID123" };

vi.mock("firebase/auth", () => {
  return {
    getAuth: () => ({
      currentUser: mockAuthUser,
    }),
  };
});

// csv parser + sha256
vi.mock("../utils/csvParser", () => ({
  parseSmilesCsv: (text: string) => {
    if (!text?.trim()) return [];
    return [
      { rowId: 1, smiles: "CCO" },
      { rowId: 2, smiles: "c1ccccc1" },
    ];
  },
}));
vi.mock("../utils/sha256", () => ({
  sha256Hex: async (_: string) => "fakehash",
}));

// Hook mocks (state derived from args)
const mockUseBatchJob = vi.fn();
const mockUseJobList = vi.fn();
const mockUseJobItems = vi.fn();

vi.mock("../hooks/useBatchJob", () => ({
  useBatchJob: (...args: any[]) => mockUseBatchJob(...args),
}));

vi.mock("../hooks/useJobList", () => ({
  useJobList: (...args: any[]) => mockUseJobList(...args),
}));

vi.mock("../hooks/useJobItems", () => ({
  useJobItems: (...args: any[]) => mockUseJobItems(...args),
}));

describe("BatchUpload (Segment 17)", () => {
  beforeEach(() => {
    vi.restoreAllMocks();
    mockUseBatchJob.mockReset();
    mockUseJobList.mockReset();
    mockUseJobItems.mockReset();

    // reset auth user each test (except the "no user" test which sets it to null)
    mockAuthUser = { uid: "testUID123" };

    // default job list: one completed job + one running job
    mockUseJobList.mockImplementation((_uid: string | null) => ({
      jobs: [
        {
          id: "job_completed",
          status: "completed",
          processed: 2,
          total: 2,
          createdAt: null,
        },
        {
          id: "job_running",
          status: "running",
          processed: 1,
          total: 5,
          createdAt: null,
        },
      ],
      loading: false,
      error: null,
    }));

    // default: no active job unless jobId present
    mockUseBatchJob.mockImplementation((_uid: string | null, jobId: string | null) => {
      if (!jobId) return { job: null, loading: false, error: null };
      return {
        job: {
          id: jobId,
          status: "completed",
          processed: 2,
          total: 2,
          successCount: 2,
          failureCount: 0,
        },
        loading: false,
        error: null,
      };
    });

    // default: items only when viewJobId present
    mockUseJobItems.mockImplementation((_uid: string | null, viewJobId: string | null) => {
      if (!viewJobId) return { items: [], loading: false, error: null };
      return {
        items: [
          {
            id: "i1",
            rowId: 1,
            smiles: "CCO",
            valid: true,
            mw: 46.07,
            logp: -0.001,
            hbd: 1,
            hba: 1,
            tpsa: 20.23,
          },
          { id: "i2", rowId: 2, smiles: "BAD", valid: false, error: "Parse error" },
        ],
        loading: false,
        error: null,
      };
    });

    // fetch mock for submit -> returns job_id
    vi.stubGlobal(
      "fetch",
      vi.fn(async () => {
        return {
          ok: true,
          json: async () => ({ job_id: "job_from_api" }),
          text: async () => "",
          status: 200,
        } as any;
      })
    );
  });

  it("renders preview after selecting a CSV file", async () => {
    render(<BatchUpload token="tkn" />);

    const input = screen.getByLabelText("csv-file") as HTMLInputElement;
    const file = new File(["smiles\nCCO\nc1ccccc1\n"], "test.csv", { type: "text/csv" });

    await userEvent.upload(input, file);

    expect(screen.getByLabelText("preview")).toBeInTheDocument();
    expect(screen.getByText("CCO")).toBeInTheDocument();
    expect(screen.getByText("c1ccccc1")).toBeInTheDocument();
  });

  it("submits, shows job_id and renders progress + View Results for completed job", async () => {
    render(<BatchUpload token="tkn" />);

    // Upload file
    const input = screen.getByLabelText("csv-file") as HTMLInputElement;
    const file = new File(["smiles\nCCO\n"], "test.csv", { type: "text/csv" });
    await userEvent.upload(input, file);

    // Submit
    await userEvent.click(screen.getByLabelText("submit"));

    // Shows job_id from API
    expect(await screen.findByLabelText("job-id")).toHaveTextContent("job_id: job_from_api");

    // Progress component should show "Completed"
    expect(screen.getByText("Completed")).toBeInTheDocument();

    // Unique completed-job button (not the list buttons)
    expect(screen.getByRole("button", { name: "View Results 📄" })).toBeInTheDocument();
  });

  it("clicking View Results shows results table and failed rows show error", async () => {
    render(<BatchUpload token="tkn" />);

    // Upload + submit to set active jobId
    const input = screen.getByLabelText("csv-file") as HTMLInputElement;
    const file = new File(["smiles\nCCO\n"], "test.csv", { type: "text/csv" });
    await userEvent.upload(input, file);
    await userEvent.click(screen.getByLabelText("submit"));

    // Click the unique completed-job View Results button
    await userEvent.click(await screen.findByRole("button", { name: "View Results 📄" }));

    // Results header appears
    expect(await screen.findByText(/Results for/i)).toBeInTheDocument();

    // Table contains a valid item + an error item
    expect(screen.getByText("Parse error")).toBeInTheDocument();
    expect(screen.getAllByText("CCO").length).toBeGreaterThanOrEqual(2);
  });

  it("job list renders past jobs, and Track switches the active job", async () => {
    render(<BatchUpload token="tkn" />);

    expect(screen.getByText("Past jobs")).toBeInTheDocument();
    expect(screen.getByText("job_completed")).toBeInTheDocument();
    expect(screen.getByText("job_running")).toBeInTheDocument();

    // Click Track for job_running -> sets jobId internally
    const trackButtons = screen.getAllByRole("button", { name: "Track" });
    await userEvent.click(trackButtons[1]); // second card = job_running

    // useBatchJob should be called with tracked job id
    expect(mockUseBatchJob).toHaveBeenCalledWith("testUID123", "job_running");
  });

  it("shows sign-in hint when uid is missing", async () => {
    mockAuthUser = null;

    render(<BatchUpload token="tkn" />);
    expect(screen.getByText(/Sign in to see job progress\/results/i)).toBeInTheDocument();
  });
});
